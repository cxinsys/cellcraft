"""
Unit tests for the plugin service layer (app/plugin/service.py) and the Docker
build helpers (app/plugin/builder.py) extracted in PR-5.

Scope: pin the extracted business logic with the filesystem and Docker/Celery
dependencies mocked (using ``tmp_path`` for real directory operations where the
rollback behavior matters). The characterization/integration tests already pin
the HTTP contract; these tests give the new service functions direct unit
coverage — especially the upload rollback path (build/DB failure must restore
files and roll back the DB).
"""
import os
import pytest
from unittest.mock import patch, MagicMock

from fastapi import HTTPException

from app.plugin import service, builder
from app.plugin.schemas import PluginCreate
from app.core.enums import PluginType


def _plugin_create(name: str = "SvcPlugin") -> PluginCreate:
    return PluginCreate(
        name=name,
        description="unit service plugin",
        author="svc_author",
        plugin_path=f"./plugin/local/{name}/",
        plugin_type=PluginType.ANALYSIS,
        dependencies={"requirements.txt": "numpy==1.24.0"},
        reference_folders=None,
        drawflow={"drawflow": {"Home": {"data": {}}}},
        rules={"0": {"name": "r", "script": "s.py"}},
        use_gpu=False,
        source="local",
        is_editable=True,
        version="1.0.0",
    )


# ---------------------------------------------------------------------------
# builder.py
# ---------------------------------------------------------------------------

class TestBuilder:
    """Docker build helpers: build-context prep and Celery dispatch."""

    def test_prepare_build_context_creates_scripts_and_dockerfile(self, tmp_path):
        """Missing scripts dir is created, a .gitkeep dummy is dropped, Dockerfile generated."""
        plugin_folder = str(tmp_path / "plug")
        os.makedirs(plugin_folder)
        script_folder = os.path.join(plugin_folder, "scripts")

        with patch("app.plugin.builder.plugin_utils.generate_plugin_dockerfile") as mock_gen:
            dockerfile_path = builder.prepare_build_context(
                plugin_folder=plugin_folder,
                script_folder=script_folder,
                use_gpu=True,
            )

        # scripts folder created and dummy .gitkeep present (was empty)
        assert os.path.isdir(script_folder)
        assert os.path.isfile(os.path.join(script_folder, ".gitkeep"))
        # Dockerfile generation delegated to utils with the expected path/flag
        assert dockerfile_path == os.path.join(plugin_folder, "Dockerfile")
        mock_gen.assert_called_once_with(plugin_folder, dockerfile_path, use_gpu=True)

    def test_prepare_build_context_skips_dummy_when_scripts_nonempty(self, tmp_path):
        """A non-empty scripts folder must not receive a .gitkeep dummy file."""
        plugin_folder = str(tmp_path / "plug")
        script_folder = os.path.join(plugin_folder, "scripts")
        os.makedirs(script_folder)
        with open(os.path.join(script_folder, "existing.py"), "w") as f:
            f.write("# real script\n")

        with patch("app.plugin.builder.plugin_utils.generate_plugin_dockerfile"):
            builder.prepare_build_context(
                plugin_folder=plugin_folder,
                script_folder=script_folder,
                use_gpu=False,
            )

        assert not os.path.exists(os.path.join(script_folder, ".gitkeep"))

    def test_dispatch_build_task_preserves_kwargs(self):
        """Celery dispatch keeps the exact wire-format kwargs the router used."""
        with patch("app.plugin.builder.build_plugin_task.apply_async") as mock_async:
            mock_async.return_value = MagicMock(id="task-xyz")
            result = builder.dispatch_build_task(plugin_name="MyPlugin", user_id=42)

        assert result.id == "task-xyz"
        mock_async.assert_called_once_with(
            args=[],
            kwargs={
                "plugin_name": "MyPlugin",
                "user_id": 42,
                "workflow_id": None,
                "algorithm_id": None,
                "task_type": "plugin_build",
            },
            ignore_result=False,
        )


# ---------------------------------------------------------------------------
# service.upload_plugin — happy path + rollback
# ---------------------------------------------------------------------------

class TestUploadPlugin:
    """upload_plugin success and rollback behavior with FS/Docker mocked."""

    @patch("app.plugin.service.invalidate_all_plugin_cache")
    @patch("app.plugin.service.crud_plugin.create_plugin")
    @patch("app.plugin.service.plugin_utils.generate_snakemake_code")
    @patch("app.plugin.service.plugin_utils.create_metadata_file")
    @patch("app.plugin.service.plugin_utils.create_dependency_folder")
    @patch("app.plugin.service.plugin_utils.create_plugin_folder")
    @patch("app.plugin.service.ensure_local_plugins_dir")
    @patch("app.plugin.service.get_plugin_path")
    @patch("app.plugin.service.os.path.exists", return_value=False)
    @patch("app.plugin.service.os.makedirs")
    def test_upload_new_plugin_success(
        self,
        _mock_makedirs,
        _mock_exists,
        mock_get_path,
        _mock_ensure,
        _mock_create_folder,
        _mock_create_dep,
        _mock_create_meta,
        _mock_gen_snake,
        mock_create_plugin,
        mock_invalidate,
    ):
        """A brand-new plugin returns the pinned success shape and commits once."""
        mock_get_path.side_effect = HTTPException(status_code=404, detail="not found")
        created = MagicMock()
        mock_create_plugin.return_value = created

        db = MagicMock()
        db.query.return_value.filter.return_value.first.return_value = None

        result = service.upload_plugin(db=db, plugin_data=_plugin_create("NewSvcPlugin"))

        assert result["message"] == "플러그인 메타데이터 업로드 성공"
        assert result["success"] is True
        assert result["redirect_to"] == "plugin_list"
        assert result["plugin"] is created
        db.commit.assert_called_once()
        mock_invalidate.assert_called_once()

    @patch("app.plugin.service.invalidate_all_plugin_cache")
    @patch("app.plugin.service.crud_plugin.create_plugin")
    @patch("app.plugin.service.plugin_utils.generate_snakemake_code")
    @patch("app.plugin.service.plugin_utils.create_metadata_file")
    @patch("app.plugin.service.plugin_utils.create_dependency_folder")
    @patch("app.plugin.service.plugin_utils.create_plugin_folder")
    @patch("app.plugin.service.ensure_local_plugins_dir")
    @patch("app.plugin.service.get_plugin_path")
    @patch("app.plugin.service.os.path.exists", return_value=False)
    @patch("app.plugin.service.os.makedirs")
    def test_upload_db_failure_rolls_back(
        self,
        _mock_makedirs,
        _mock_exists,
        mock_get_path,
        _mock_ensure,
        _mock_create_folder,
        _mock_create_dep,
        _mock_create_meta,
        _mock_gen_snake,
        mock_create_plugin,
        _mock_invalidate,
    ):
        """A DB write failure rolls the transaction back and yields HTTPException(500)."""
        mock_get_path.side_effect = HTTPException(status_code=404, detail="not found")
        mock_create_plugin.side_effect = Exception("DB write boom")

        db = MagicMock()
        db.query.return_value.filter.return_value.first.return_value = None

        with pytest.raises(HTTPException) as exc:
            service.upload_plugin(db=db, plugin_data=_plugin_create("RollbackSvcPlugin"))

        assert exc.value.status_code == 500
        assert "플러그인 메타데이터 업데이트 실패" in exc.value.detail
        # DB was rolled back, never committed.
        db.rollback.assert_called_once()
        db.commit.assert_not_called()

    @patch("app.plugin.service.ensure_local_plugins_dir")
    @patch("app.plugin.service.get_plugin_path")
    def test_upload_over_official_plugin_is_bypassed(self, mock_get_path, _mock_ensure):
        """ANOMALY (pinned): the official-plugin 403 guard is silently swallowed.

        get_plugin_path resolves to an official plugin, so the guard raises
        HTTPException(403); but upload_plugin catches HTTPException and passes,
        so the upload proceeds. Here it fails later at folder creation and
        surfaces a 500 — matching the integration characterization test.
        """
        mock_get_path.return_value = ("./plugin/official/Off/", "official")

        db = MagicMock()
        db.query.return_value.filter.return_value.first.return_value = None

        with patch("app.plugin.service.plugin_utils.create_plugin_folder",
                   side_effect=Exception("boom")):
            with pytest.raises(HTTPException) as exc:
                service.upload_plugin(db=db, plugin_data=_plugin_create("Off"))

        # Not a 403 — the guard was bypassed and a later failure produced 500.
        assert exc.value.status_code == 500


class TestRollbackUpload:
    """Direct coverage of the _rollback_upload step function (real filesystem)."""

    def test_rollback_restores_backup(self, tmp_path):
        """When a backup exists, the plugin folder is removed and restored from backup."""
        plugin_folder = str(tmp_path / "plugin")
        backup_folder = str(tmp_path / "backup")
        os.makedirs(plugin_folder)
        os.makedirs(backup_folder)
        # A file only present in the backup proves restoration happened.
        with open(os.path.join(backup_folder, "orig.txt"), "w") as f:
            f.write("original")
        # A file only in the (failed) new folder proves it was discarded.
        with open(os.path.join(plugin_folder, "new.txt"), "w") as f:
            f.write("new")

        service._rollback_upload(plugin_folder=plugin_folder, backup_folder=backup_folder)

        assert os.path.isfile(os.path.join(plugin_folder, "orig.txt"))
        assert not os.path.exists(os.path.join(plugin_folder, "new.txt"))
        # Backup folder is consumed by the restore.
        assert not os.path.exists(backup_folder)

    def test_rollback_without_backup_removes_new_folder(self, tmp_path):
        """With no backup, a fresh-install failure simply removes the new folder."""
        plugin_folder = str(tmp_path / "plugin")
        os.makedirs(plugin_folder)
        with open(os.path.join(plugin_folder, "new.txt"), "w") as f:
            f.write("new")

        service._rollback_upload(plugin_folder=plugin_folder, backup_folder=None)

        assert not os.path.exists(plugin_folder)


# ---------------------------------------------------------------------------
# CRUD-delegating services
# ---------------------------------------------------------------------------

class TestSimpleServices:
    """Thin CRUD-delegating service functions."""

    @patch("app.plugin.service.crud_plugin.get_plugin_by_name")
    def test_get_plugin_info_found(self, mock_get):
        plugin = MagicMock()
        mock_get.return_value = plugin
        db = MagicMock()

        result = service.get_plugin_info(db=db, plugin_name="P")

        assert result == {"plugin": plugin}

    @patch("app.plugin.service.crud_plugin.get_plugin_by_name")
    def test_get_plugin_info_not_found_maps_to_500(self, mock_get):
        """ANOMALY (pinned): the inner 404 is re-wrapped as 500 by the broad except.

        get_plugin_info raises HTTPException(404) when the plugin is missing, but
        the surrounding ``except Exception`` catches it and re-raises as 500 —
        matching the current router behavior.
        """
        mock_get.return_value = None
        db = MagicMock()

        with pytest.raises(HTTPException) as exc:
            service.get_plugin_info(db=db, plugin_name="Missing")

        assert exc.value.status_code == 500

    @patch("app.plugin.service.invalidate_all_plugin_cache")
    @patch("app.plugin.service.crud_plugin.get_plugin_by_id")
    @patch("app.plugin.service.crud_plugin.associate_user_plugin")
    def test_associate_plugin(self, mock_assoc, mock_get_by_id, mock_invalidate):
        plugin = MagicMock()
        mock_get_by_id.return_value = plugin
        db = MagicMock()
        user = MagicMock(id=7)

        result = service.associate_plugin(db=db, user=user, plugin_id=3)

        mock_assoc.assert_called_once_with(db, 7, 3)
        assert result == {"message": "Plugin associated with user successfully", "plugin": plugin}
        mock_invalidate.assert_called_once()

    @patch("app.plugin.service.invalidate_all_plugin_cache")
    @patch("app.plugin.service.crud_plugin.get_plugin_by_id")
    @patch("app.plugin.service.crud_plugin.dissociate_user_plugin")
    def test_dissociate_plugin(self, mock_dissoc, mock_get_by_id, mock_invalidate):
        plugin = MagicMock()
        mock_get_by_id.return_value = plugin
        db = MagicMock()
        user = MagicMock(id=9)

        result = service.dissociate_plugin(db=db, user=user, plugin_id=5)

        mock_dissoc.assert_called_once_with(db, 9, 5)
        assert result == {"message": "Plugin dissociated from user successfully", "plugin": plugin}
        mock_invalidate.assert_called_once()


class TestBuildStatusService:
    """Async build status/cancel services (Celery mocked)."""

    @pytest.mark.asyncio
    @patch("app.plugin.service.AsyncResult")
    async def test_get_build_status_pending(self, mock_async_result):
        mock_async_result.return_value = MagicMock(state="PENDING", info=None)

        result = await service.get_build_status(task_id="t1")

        assert result["task_id"] == "t1"
        assert result["state"] == "PENDING"
        assert result["info"] == {"message": "Build is waiting to start..."}

    @pytest.mark.asyncio
    @patch("app.plugin.service.invalidate_all_plugin_cache")
    @patch("app.plugin.service.AsyncResult")
    @patch("app.plugin.service.celery_app")
    async def test_cancel_build_task(self, mock_celery, mock_async_result, _mock_invalidate):
        mock_async_result.return_value = MagicMock(state="REVOKED")

        result = await service.cancel_build_task(task_id="t2")

        mock_celery.control.revoke.assert_called_once_with("t2", terminate=True)
        assert result["task_id"] == "t2"
        assert result["message"] == "Build task t2 cancellation requested"
        assert result["state"] == "REVOKED"
