"""
Characterization tests for the file upload + download-security surface.

Purpose: pin the CURRENT behavior of the files endpoints so the later service +
storage/security extraction refactor (PR-8) can prove behavior is unchanged.

Endpoints covered (app/file/router.py):
- POST /routes/files/upload           (fileUpload; lines 25-132)
- GET  /routes/files/result/{filename} (download_result_file; lines 538-555)
    -> exercises file_security.validate_file_path path-traversal rejection (403)

Source of truth read while writing these tests:
- app/file/router.py :: fileUpload, download_result_file
- app/file/security.py :: validate_file_upload, validate_file_path
- app/core/config.py :: ALLOWED_EXTENSIONS ({.h5ad,.csv,.json,.txt}),
  MAX_FILES_PER_UPLOAD (20)

Pinned behavior:
- Successful single upload returns the created File record (raw model fields:
  id, file_name, file_size, file_path, folder, user_id) and a DB row exists.
  Filename format "folder_actualname.ext" is split on the FIRST underscore.
- Disallowed extension -> 400 "Invalid file type".
- Empty file -> 400 "Empty file not allowed".
- Duplicate filename (already in DB) -> 400 "already exists".
- Path traversal on /result/{filename} -> 403 "Path traversal attempt detected".
"""
import os
import shutil

import pytest
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app import models
from app.file import crud as crud_file
from app.core.config import settings

UPLOAD_URL = f"{settings.ROUTES_STR}/files/upload"
RESULT_URL = f"{settings.ROUTES_STR}/files/result"


@pytest.mark.integration
@pytest.mark.characterization
class TestCharacterizationFileUpload:
    """Freeze current upload success/reject behavior and DB persistence."""

    def test_upload_single_file_success_record_and_db(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict,
        sample_user: models.User,
        db_session: Session,
    ):
        """A valid .h5ad upload returns the raw File record and persists a DB row."""
        upload_file = upload_file_factory(
            filename="charfolder_char_data.h5ad",
            content=b"mock_h5ad_content" * 50,
        )

        try:
            response = client.post(
                UPLOAD_URL,
                headers=auth_headers,
                files={"files": (
                    upload_file.filename,
                    upload_file.file,
                    upload_file.content_type,
                )},
            )

            assert response.status_code == 200, response.text
            data = response.json()

            # Endpoint returns the raw SQLAlchemy File object for a single upload.
            # Filename splits on the FIRST underscore: folder="charfolder",
            # file_name="char_data.h5ad".
            assert data["file_name"] == "char_data.h5ad"
            assert data["folder"] == "charfolder"
            assert data["user_id"] == sample_user.id
            assert "id" in data and "file_size" in data and "file_path" in data

            # DB row created.
            db_file = crud_file.get_user_file(db_session, sample_user.id, "char_data.h5ad")
            assert db_file is not None
            assert db_file.folder == "charfolder"
        finally:
            # The endpoint writes to the real relative path ./user/{username}/data.
            user_dir = f"./user/{sample_user.username}"
            if os.path.exists(user_dir):
                shutil.rmtree(user_dir, ignore_errors=True)

    def test_upload_disallowed_extension_rejected_400(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
    ):
        """A .sh file is rejected by the extension whitelist with 400."""
        upload_file = upload_file_factory(
            filename="charfolder_evil.sh",
            content=b"#!/bin/sh\necho hi\n" * 5,
            content_type="application/x-sh",
        )

        response = client.post(
            UPLOAD_URL,
            headers=auth_headers,
            files={"files": (
                upload_file.filename,
                upload_file.file,
                upload_file.content_type,
            )},
        )

        assert response.status_code == 400
        assert "invalid file type" in response.json()["detail"].lower()

    def test_upload_empty_file_rejected_400(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
    ):
        """A zero-byte file is rejected with 400 'Empty file not allowed'."""
        upload_file = upload_file_factory(
            filename="charfolder_empty.h5ad",
            content=b"",
        )

        response = client.post(
            UPLOAD_URL,
            headers=auth_headers,
            files={"files": (
                upload_file.filename,
                upload_file.file,
                upload_file.content_type,
            )},
        )

        assert response.status_code == 400
        assert "empty file" in response.json()["detail"].lower()

    def test_upload_duplicate_filename_rejected_400(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        sample_file_metadata: models.File,
    ):
        """Uploading a filename already present in the DB is rejected with 400."""
        # sample_file_metadata.file_name is "test_data.h5ad"; prefix a folder so the
        # split yields the same actual filename.
        upload_file = upload_file_factory(
            filename=f"charfolder_{sample_file_metadata.file_name}",
            content=b"new_content" * 50,
        )

        response = client.post(
            UPLOAD_URL,
            headers=auth_headers,
            files={"files": (
                upload_file.filename,
                upload_file.file,
                upload_file.content_type,
            )},
        )

        assert response.status_code == 400
        assert "already exists" in response.json()["detail"]

    def test_upload_requires_authentication(
        self,
        client: TestClient,
        upload_file_factory,
    ):
        """Upload without auth headers is rejected with 401."""
        upload_file = upload_file_factory(
            filename="charfolder_x.h5ad",
            content=b"data" * 10,
        )
        response = client.post(
            UPLOAD_URL,
            files={"files": (
                upload_file.filename,
                upload_file.file,
                upload_file.content_type,
            )},
        )
        assert response.status_code == 401


@pytest.mark.integration
@pytest.mark.characterization
class TestCharacterizationFileSecurityRejection:
    """Freeze file_security path-traversal rejection on the result download route."""

    def test_download_result_path_traversal_returns_404(
        self,
        client: TestClient,
        auth_headers: dict,
    ):
        """URL-encoded path traversal on /result returns 404 today.

        ANOMALY (pinned): Starlette normalises the URL-encoded ``%2E%2E%2F``
        segments before the path parameter is matched, so the traversal string
        never reaches ``validate_file_path``. The router returns 404 (no route
        match) instead of the intended 403. The security check is therefore
        unreachable via URL-encoded traversal at this endpoint today.
        """
        response = client.get(
            f"{RESULT_URL}/%2E%2E%2F%2E%2E%2Fetc%2Fpasswd",
            headers=auth_headers,
        )

        assert response.status_code == 404
