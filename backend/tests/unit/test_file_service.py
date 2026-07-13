"""
Unit tests for the file service layer (app/file/service.py) extracted in PR-8.

Scope: pin the extracted business logic with the filesystem / crud / security
boundaries mocked. The characterization + integration tests already pin the HTTP
contract; these give the new service functions direct coverage and assert the
domain-exception mapping (``ValidationFailedError`` -> 400, ``NotFoundError`` ->
404) that the global handler renders as the unchanged ``{"detail": ...}`` body.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.core.exceptions import NotFoundError, ValidationFailedError
from app.file import service


def _user(username="filesvc", user_id=3):
    u = MagicMock()
    u.id = user_id
    u.username = username
    return u


def _file_info(file_name="d.h5ad", source="local", anno_column="cell_type",
               option_file_name="opt"):
    fi = MagicMock()
    fi.file_name = file_name
    fi.source = source
    fi.anno_column = anno_column
    fi.option_file_name = option_file_name
    return fi


# ---------------------------------------------------------------------------
# save_upload: validation branches (async)
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
class TestSaveUpload:
    async def test_too_many_files_raises_400(self):
        db = MagicMock()
        files = [MagicMock(), MagicMock(), MagicMock()]
        with patch("app.core.config.MAX_FILES_PER_UPLOAD", 1):
            with pytest.raises(ValidationFailedError) as exc:
                await service.save_upload(db=db, files=files, current_user=_user())
        assert exc.value.status_code == 400
        assert "Too many files" in exc.value.detail

    async def test_invalid_filename_format_raises_400(self):
        db = MagicMock()
        files = [MagicMock()]
        with patch("app.core.config.MAX_FILES_PER_UPLOAD", 20), \
                patch("app.file.service.os.makedirs"), \
                patch("app.file.security.validate_file_upload", return_value="nounderscore.h5ad"):
            with pytest.raises(ValidationFailedError) as exc:
                await service.save_upload(db=db, files=files, current_user=_user())
        assert exc.value.status_code == 400
        assert "Invalid filename format" in exc.value.detail

    async def test_duplicate_filename_raises_400(self):
        db = MagicMock()
        files = [MagicMock()]
        with patch("app.core.config.MAX_FILES_PER_UPLOAD", 20), \
                patch("app.file.service.os.makedirs"), \
                patch("app.file.security.validate_file_upload", return_value="folder_dup.h5ad"), \
                patch("app.file.service.crud_file.get_user_file", return_value=MagicMock()):
            with pytest.raises(ValidationFailedError) as exc:
                await service.save_upload(db=db, files=files, current_user=_user())
        assert exc.value.status_code == 400
        assert "already exists" in exc.value.detail


# ---------------------------------------------------------------------------
# Simple lookup functions
# ---------------------------------------------------------------------------

class TestLookups:
    def test_find_user_file_found(self):
        db = MagicMock()
        rec = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=rec):
            assert service.find_user_file(db=db, current_user=_user(), file_info=_file_info()) is rec

    def test_find_user_file_missing_raises_400(self):
        db = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=None):
            with pytest.raises(ValidationFailedError) as exc:
                service.find_user_file(db=db, current_user=_user(), file_info=_file_info())
        assert exc.value.status_code == 400
        assert exc.value.detail == "this file not exists in your files"

    def test_find_user_folder_missing_raises_400(self):
        db = MagicMock()
        folder = MagicMock()
        folder.folder_name = "nope"
        with patch("app.file.service.crud_file.get_user_folder", return_value=None):
            with pytest.raises(ValidationFailedError) as exc:
                service.find_user_folder(db=db, current_user=_user(), folder=folder)
        assert exc.value.detail == "this folder not exists in your folders"

    def test_update_file_missing_raises_400(self):
        db = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=None):
            with pytest.raises(ValidationFailedError):
                service.update_file(db=db, current_user=_user(), file_info=_file_info())

    def test_read_user_folder_shapes_walk(self):
        walk = [("./user/filesvc/data", [], ["a.h5ad", "b.csv"])]
        with patch("app.file.service.os.walk", return_value=walk):
            res = service.read_user_folder(current_user=_user())
        assert res == [("data", ["a.h5ad", "b.csv"])]


# ---------------------------------------------------------------------------
# delete_file
# ---------------------------------------------------------------------------

class TestDeleteFile:
    def test_missing_raises_400(self):
        db = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=None):
            with pytest.raises(ValidationFailedError) as exc:
                service.delete_file(db=db, current_user=_user(), file_info=_file_info())
        assert exc.value.detail == "this file not exists in your files"

    def test_deletes_and_returns_crud_result(self):
        db = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=MagicMock()), \
                patch("app.file.security.validate_file_path"), \
                patch("app.file.service.os.path.exists", return_value=False), \
                patch("app.file.service.crud_file.delete_user_file", return_value="deleted"):
            out = service.delete_file(db=db, current_user=_user(), file_info=_file_info())
        assert out == "deleted"


# ---------------------------------------------------------------------------
# h5ad columns/clusters: ownership + not-found + cache
# ---------------------------------------------------------------------------

class TestH5adColumns:
    def test_unowned_local_file_raises_400(self):
        db = MagicMock()
        with patch("app.file.service.crud_file.get_user_file", return_value=None):
            with pytest.raises(ValidationFailedError):
                service.get_h5ad_columns(db=db, current_user=_user(), file_info=_file_info(source="local"))

    def test_missing_file_raises_404(self):
        db = MagicMock()
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/x/y.h5ad"), \
                patch("app.file.service.os.path.isfile", return_value=False):
            with pytest.raises(NotFoundError) as exc:
                service.get_h5ad_columns(db=db, current_user=_user(), file_info=_file_info(source="shared"))
        assert exc.value.status_code == 404
        assert exc.value.detail == "File not found"

    def test_cache_hit_returns_cached(self):
        db = MagicMock()
        cached = {"anno_columns": ["a"], "pseudo_columns": []}
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/x/y.h5ad"), \
                patch("app.file.service.os.path.isfile", return_value=True), \
                patch("app.file.service.get_cached_columns", return_value=cached):
            out = service.get_h5ad_columns(db=db, current_user=_user(), file_info=_file_info(source="shared"))
        assert out == cached

    def test_clusters_cache_hit(self):
        db = MagicMock()
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/x/y.h5ad"), \
                patch("app.file.service.os.path.isfile", return_value=True), \
                patch("app.file.service.get_cached_clusters", return_value=["c1", "c2"]):
            out = service.get_h5ad_clusters(db=db, current_user=_user(), file_info=_file_info(source="shared"))
        assert out == {"clusters": ["c1", "c2"]}


# ---------------------------------------------------------------------------
# setup / result / shared / download helpers
# ---------------------------------------------------------------------------

class TestFileAccess:
    def test_read_setup_missing_raises_404(self):
        with patch("app.file.security.validate_file_path"), \
                patch("app.file.service.os.path.exists", return_value=False):
            with pytest.raises(NotFoundError) as exc:
                service.read_setup_file(current_user=_user(), file_name="x.json")
        assert exc.value.status_code == 404

    def test_check_convert_missing_raises_400(self):
        with patch("app.file.service.os.path.isfile", return_value=False):
            with pytest.raises(ValidationFailedError):
                service.check_convert(current_user=_user(), file_name="d.h5ad")

    def test_download_result_missing_raises_404(self):
        with patch("app.file.security.validate_file_path"), \
                patch("app.file.service.os.path.exists", return_value=False):
            with pytest.raises(NotFoundError):
                service.download_result_file(current_user=_user(), filename="r.csv")

    def test_download_tutorial_missing_raises_404(self):
        with patch("app.file.security.validate_file_path"), \
                patch("app.file.service.os.path.exists", return_value=False):
            with pytest.raises(NotFoundError):
                service.download_tutorial_file(current_user=_user(), filename="t.html")

    def test_get_shared_files_no_dir_returns_empty(self):
        with patch("app.file.service.os.path.isdir", return_value=False):
            assert service.get_shared_files() == []

    def test_find_shared_missing_raises_404(self):
        with patch("app.file.security.validate_file_path"), \
                patch("app.file.service.os.path.exists", return_value=False):
            with pytest.raises(NotFoundError) as exc:
                service.find_shared_file(file_info=_file_info())
        assert exc.value.detail == "Shared file not found"


@pytest.mark.asyncio
class TestDownloadData:
    async def test_missing_raises_400(self):
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/x/y.csv"), \
                patch("app.file.service.os.path.isfile", return_value=False):
            with pytest.raises(ValidationFailedError) as exc:
                await service.download_data_file(current_user=_user(), filename="y.csv", source=None)
        assert exc.value.status_code == 400
        assert exc.value.detail == "this file not exists in your files"
