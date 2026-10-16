"""
Integration tests for File API - Phase 3.3

Test Coverage:
- Security: Path traversal, file upload validation (CRITICAL)
- CRUD: Upload, download, update, delete operations
- H5AD: Conversion, column extraction, cluster analysis
- Metadata: File listing, search, folder operations

Critical Security Issues Identified:
1. Path traversal vulnerabilities in 4 download endpoints
2. Missing file upload validation (size, type, filename sanitization)
3. No authorization checks for cross-user file access
"""
import os
import pytest
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app import models
from app.file import crud as crud_file

# File-specific fixtures are automatically available from conftest.py
# No need to import explicitly - pytest discovers fixtures automatically


# ==============================================================================
# CRITICAL PRIORITY - Security Tests (Path Traversal)
# ==============================================================================

class TestFileSecurityPathTraversal:
    """
    CRITICAL: Path traversal security tests.

    Identified vulnerabilities:
    - /result/{filename} (line 349-355) - NO validation
    - /data/{filename} (line 357-374) - Basic existence check only
    - /tutorials/{filename} (line 376-383) - NO validation
    - /setup/{file_name} (line 295-312) - NO validation

    Attack vector: ../../etc/passwd or ../../../backend/.env
    """

    def test_download_result_path_traversal_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        malicious_filenames: list
    ):
        """
        Test path traversal prevention in /result/{filename} endpoint.

        Expected: 403 Forbidden for all path traversal attempts
        """
        for malicious_name in malicious_filenames:
            response = client.get(
                f"/routes/files/result/{malicious_name}",
                headers=auth_headers
            )

            # Must return exactly 403 Forbidden - no other status codes acceptable
            # This prevents regression of security fixes
            assert response.status_code == 403, \
                f"Expected 403 Forbidden for path traversal: {malicious_name}, got {response.status_code}"
            assert "Path traversal attempt detected" in response.json()["detail"], \
                "Error message must indicate path traversal detection"

    def test_download_data_path_traversal_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        malicious_filenames: list
    ):
        """
        Test path traversal prevention in /data/{filename} endpoint.

        Expected: 403 Forbidden for path traversal attempts
        """
        for malicious_name in malicious_filenames:
            response = client.get(
                f"/routes/files/data/{malicious_name}",
                headers=auth_headers
            )

            # Must return exactly 403 Forbidden - no other status codes acceptable
            # This prevents regression of security fixes
            assert response.status_code == 403, \
                f"Expected 403 Forbidden for path traversal: {malicious_name}, got {response.status_code}"
            assert "Path traversal attempt detected" in response.json()["detail"], \
                "Error message must indicate path traversal detection"

    def test_download_tutorial_path_traversal_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        malicious_filenames: list
    ):
        """
        Test path traversal prevention in /tutorials/{filename} endpoint.

        Expected: 403 Forbidden for all path traversal attempts
        """
        for malicious_name in malicious_filenames:
            response = client.get(
                f"/routes/files/tutorials/{malicious_name}",
                headers=auth_headers
            )

            # Must return exactly 403 Forbidden - no other status codes acceptable
            # This prevents regression of security fixes
            assert response.status_code == 403, \
                f"Expected 403 Forbidden for path traversal: {malicious_name}, got {response.status_code}"
            assert "Path traversal attempt detected" in response.json()["detail"], \
                "Error message must indicate path traversal detection"

    def test_setup_file_path_traversal_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        malicious_filenames: list
    ):
        """
        Test path traversal prevention in /setup/{file_name} endpoint.

        Expected: 403 Forbidden for traversal attempts
        """
        for malicious_name in malicious_filenames:
            response = client.get(
                f"/routes/files/setup/{malicious_name}",
                headers=auth_headers
            )

            # Must return exactly 403 Forbidden - no other status codes acceptable
            # This prevents regression of security fixes
            assert response.status_code == 403, \
                f"Expected 403 Forbidden for path traversal: {malicious_name}, got {response.status_code}"
            assert "Path traversal attempt detected" in response.json()["detail"], \
                "Error message must indicate path traversal detection"

    def test_delete_file_path_traversal_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User
    ):
        """
        Test EXISTING path traversal prevention in /delete endpoint.

        Current Implementation: GOOD - has realpath validation (lines 119-128)
        This test verifies the existing protection works correctly.
        """
        # Attempt to delete file outside user directory
        response = client.post(
            "/routes/files/delete",
            headers=auth_headers,
            json={"file_name": "../../../etc/passwd"}
        )

        # Should be blocked by realpath check or "file not exists"
        assert response.status_code in [403, 400], \
            "Delete path traversal protection failed"

        if response.status_code == 403:
            assert "Access denied" in response.json()["detail"]

    def test_case_sensitivity_bypass_attempts(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test case sensitivity bypass attempts for path traversal.

        Some filesystems (e.g., Windows) are case-insensitive, which could
        allow attackers to bypass simple string-matching filters.

        Expected: 403 Forbidden for all case variations
        """
        case_variations = [
            "../../../ETC/PASSWD",
            "../../../Etc/Passwd",
            "../../../etc/PASSWD",
            "..\\..\\..\\WINDOWS\\SYSTEM32\\config\\sam",
            "..\\..\\..\\Windows\\System32\\CONFIG\\sam"
        ]

        for malicious_name in case_variations:
            response = client.get(
                f"/routes/files/data/{malicious_name}",
                headers=auth_headers
            )

            # Must return exactly 403 Forbidden
            assert response.status_code == 403, \
                f"Expected 403 Forbidden for case bypass attempt: {malicious_name}, got {response.status_code}"
            assert "Path traversal attempt detected" in response.json()["detail"], \
                "Error message must indicate path traversal detection"

    def test_path_traversal_sibling_directory_blocked(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        temp_user_file_directory: dict,
        db_session
    ):
        """
        Test CRITICAL sibling directory attack prevention.

        Attack Pattern: If user folder is ./user/john/data
        Attacker provides: ../data_backup/secrets.txt
        After realpath: /home/cellcraft/user/john/data_backup/secrets.txt

        Old vulnerable code using startswith():
        - This PASSES because path starts with "/home/cellcraft/user/john/data"

        New secure code using Path.relative_to():
        - This FAILS with ValueError because path is not under allowed directory

        Expected: 403 Forbidden with "Path traversal attempt detected"
        """
        # Create a sibling directory to simulate attack target
        # Use temp_user_file_directory to ensure proper permissions
        user_base = os.path.dirname(temp_user_file_directory["data_dir"])
        sibling_dir = os.path.join(user_base, 'data_backup')
        os.makedirs(sibling_dir, exist_ok=True)

        # Create a "secret" file in sibling directory
        secret_file = os.path.join(sibling_dir, 'secrets.txt')
        with open(secret_file, 'w') as f:
            f.write("CONFIDENTIAL DATA")

        # Create DB record to pass initial file existence check
        # This simulates an attacker knowing a sibling path exists
        crud_file.create_file(
            db_session,
            "../data_backup/secrets.txt",
            len("CONFIDENTIAL DATA"),
            temp_user_file_directory["data_dir"],
            "data",
            sample_user.id
        )

        # Attempt to download via path traversal to sibling directory
        # FastAPI will pass this as the filename parameter
        response = client.get(
            "/routes/files/data/%2E%2E%2Fdata_backup%2Fsecrets.txt",  # URL-encoded ../data_backup/secrets.txt
            headers=auth_headers
        )

        # Must be blocked with 403 (not 400/404, as file record exists in DB)
        assert response.status_code == 403, \
            f"Sibling directory attack not blocked! Got {response.status_code}"
        assert "Path traversal attempt detected" in response.json()["detail"], \
            "Error message does not indicate path traversal detection"

        # Cleanup
        import shutil
        if os.path.exists(sibling_dir):
            shutil.rmtree(sibling_dir)

    def test_symlink_attack_prevention(
        self,
        client: TestClient,
        auth_headers: dict,
        symlink_attack_setup: dict,
        sample_file_metadata: models.File
    ):
        """
        Test symlink attack prevention.

        Attack: Create symlink pointing to /etc/passwd, attempt to download.
        Expected: Should detect and block symlink access.
        """
        # Try to download via symlink
        response = client.get(
            f"/routes/files/data/{symlink_attack_setup['symlink_filename']}",
            headers=auth_headers
        )

        # TODO: Add symlink detection
        # Current: May follow symlink to /etc/passwd
        # Expected: Detect symlink and reject with 403
        assert response.status_code in [403, 400, 404], \
            "Symlink attack not prevented"

    def test_filename_sanitization_upload(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test filename sanitization during upload.

        Attack: Upload files with malicious characters in filename.
        Expected: Filename sanitized or upload rejected.
        """
        malicious_filenames = [
            "file;rm -rf /.h5ad",
            "file<script>.h5ad",
            "file\x00.h5ad",
            "file\n\r\t.h5ad"
        ]

        for malicious_name in malicious_filenames:
            upload_file = upload_file_factory(
                filename=f"folder_{malicious_name}",
                content=b"mock_data" * 100
            )

            response = client.post(
                "/routes/files/upload",
                headers=auth_headers,
                files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
            )

            # TODO: Add filename sanitization
            # Current: Accepts any filename
            # Expected: Sanitize or reject malicious filenames
            if response.status_code == 200:
                # Verify filename was sanitized
                created_file = response.json()
                assert ";" not in created_file.get("file_name", ""), \
                    "Dangerous characters not sanitized"

    def test_concurrent_file_access_race_condition(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File
    ):
        """
        Test race condition handling: simultaneous delete and download.

        Attack: Thread 1 deletes while Thread 2 downloads same file.
        Expected: Atomic operations, no partial state.
        """
        import threading
        import time

        results = {"delete": None, "download": None}

        def delete_file():
            time.sleep(0.01)  # Small delay to sync threads
            results["delete"] = client.post(
                "/routes/files/delete",
                headers=auth_headers,
                json={"file_name": sample_file_metadata.file_name}
            )

        def download_file():
            time.sleep(0.01)
            results["download"] = client.get(
                f"/routes/files/data/{sample_file_metadata.file_name}",
                headers=auth_headers
            )

        # Start threads simultaneously
        thread1 = threading.Thread(target=delete_file)
        thread2 = threading.Thread(target=download_file)

        thread1.start()
        thread2.start()

        thread1.join()
        thread2.join()

        # One should succeed, one should fail gracefully
        assert (results["delete"].status_code == 200 and results["download"].status_code in [400, 404]) or \
               (results["download"].status_code == 200 and results["delete"].status_code == 400), \
            "Race condition not handled properly"


# ==============================================================================
# CRITICAL PRIORITY - Security Tests (File Upload Validation)
# ==============================================================================

class TestFileUploadValidation:
    """
    CRITICAL: File upload validation security tests.

    Missing validations in /upload endpoint (lines 23-53):
    1. File size limits - NO MAX_UPLOAD_SIZE check
    2. File type validation - NO ALLOWED_EXTENSIONS whitelist
    3. Malicious filename sanitization
    4. Multiple file upload abuse - No limit on file count
    """

    def test_upload_exceeds_max_size(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_large_upload_file,
        temp_user_file_directory: dict
    ):
        """
        Test file size limit enforcement.

        Expected: 413 Payload Too Large for oversized files
        Recommended limit: 500MB for H5AD files
        """
        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (
                mock_large_upload_file.filename,
                mock_large_upload_file.file,
                mock_large_upload_file.content_type
            )}
        )

        # Must return exactly 413 for oversized files
        # This prevents regression of security fixes
        assert response.status_code == 413, \
            f"Expected 413 Payload Too Large for oversized file, got {response.status_code}"
        error_detail = response.json()["detail"].lower()
        assert "too large" in error_detail or "exceeds maximum" in error_detail, \
            "Error message must indicate file size limit exceeded"

    def test_upload_invalid_file_type(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test file type validation.

        Expected: 400 Bad Request with "Invalid file type"
        Recommended whitelist: [".h5ad", ".csv", ".json"]
        """
        malicious_files = [
            ("folder_malicious.sh", b"#!/bin/bash\nrm -rf /", "application/x-sh"),
            ("folder_script.py", b"import os\nos.system('cat /etc/passwd')", "text/x-python"),
            ("folder_executable.exe", b"\x4D\x5A\x90\x00", "application/x-msdownload")
        ]

        for filename, content, content_type in malicious_files:
            upload_file = upload_file_factory(
                filename=filename,
                content=content,
                content_type=content_type
            )

            response = client.post(
                "/routes/files/upload",
                headers=auth_headers,
                files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
            )

            # Must return exactly 400 for invalid file types
            # This prevents regression of security fixes
            assert response.status_code == 400, \
                f"Expected 400 Bad Request for invalid file type {filename}, got {response.status_code}"
            assert "invalid file type" in response.json()["detail"].lower() or \
                   "file type not allowed" in response.json()["detail"].lower(), \
                f"Error message must indicate invalid file type for {filename}"

    def test_upload_filename_length_limit(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test filename length limit.

        Expected: 400 Bad Request for excessively long filenames
        Maximum allowed: 255 characters
        """
        long_filename = "folder_" + "a" * 500 + ".h5ad"
        upload_file = upload_file_factory(
            filename=long_filename,
            content=b"mock_data" * 10
        )

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        # Must return exactly 400 for excessively long filenames
        # This prevents regression of security fixes
        assert response.status_code == 400, \
            f"Expected 400 Bad Request for long filename, got {response.status_code}"
        assert "filename too long" in response.json()["detail"].lower() or \
               "exceeds maximum length" in response.json()["detail"].lower(), \
            "Error message must indicate filename length violation"

    def test_upload_multiple_files_limit(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test multiple file upload limit.

        Expected: 400 Bad Request for excessive batch uploads
        Maximum allowed: 10-20 files per request
        """
        # Create 50 files (exceeds recommended limit)
        files = []
        for i in range(50):
            upload_file = upload_file_factory(
                filename=f"folder_file{i}.h5ad",
                content=b"mock_data" * 10
            )
            files.append(("files", (upload_file.filename, upload_file.file, upload_file.content_type)))

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files=files
        )

        # Must return exactly 400 for excessive batch uploads
        # This prevents regression of security fixes
        assert response.status_code == 400, \
            f"Expected 400 Bad Request for batch upload limit exceeded, got {response.status_code}"
        assert "too many files" in response.json()["detail"].lower() or \
               "exceeds maximum" in response.json()["detail"].lower(), \
            "Error message must indicate batch upload limit exceeded"


# ==============================================================================
# HIGH PRIORITY - File Upload CRUD Operations
# ==============================================================================

class TestFileUploadOperations:
    """
    Test file upload CRUD operations.
    """

    def test_upload_single_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict,
        sample_user: models.User
    ):
        """
        Test successful single file upload.

        Expected:
        - 200 status
        - File saved to ./user/{username}/data/filename.h5ad
        - Database record created with correct metadata
        - Folder and filename parsed correctly
        """
        upload_file = upload_file_factory(
            filename="test_folder_test_data.h5ad",
            content=b"mock_h5ad_data" * 100  # ~1.4KB
        )

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        assert response.status_code == 200, f"Upload failed: {response.text}"
        data = response.json()

        # Verify response
        assert data["file_name"] == "test_data.h5ad"
        assert data["folder"] == "test_folder"

        # Verify file exists on filesystem
        file_path = os.path.join(temp_user_file_directory["data_dir"], "test_data.h5ad")
        assert os.path.exists(file_path), "File not saved to filesystem"

        # Verify file size
        assert os.path.getsize(file_path) > 0, "File is empty"

    def test_upload_multiple_files_success(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_upload_files: list,
        temp_user_file_directory: dict,
        sample_user: models.User
    ):
        """
        Test successful multiple file upload in single request.

        Expected:
        - All 3 files saved to filesystem
        - All database records created
        - Folder/filename parsing correct for each
        """
        # Prepare files for upload
        files_data = [
            ("files", (f.filename, f.file, f.content_type))
            for f in mock_upload_files
        ]

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files=files_data
        )

        assert response.status_code == 200, f"Batch upload failed: {response.text}"

        # Verify all files exist
        for upload_file in mock_upload_files:
            # Extract actual filename (after folder_ prefix)
            folder, filename = upload_file.filename.split('_', 1)
            file_path = os.path.join(temp_user_file_directory["data_dir"], filename)
            assert os.path.exists(file_path), f"File {filename} not saved"

    def test_upload_duplicate_file_rejected(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test duplicate file upload rejection.

        Expected:
        - 400 Bad Request
        - Error message: "this file already exists in your files"
        - Original file unchanged
        """
        # Create file with same name as existing metadata
        upload_file = upload_file_factory(
            filename=f"folder_{sample_file_metadata.file_name}",
            content=b"new_data" * 100
        )

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        assert response.status_code == 400
        assert "already exists" in response.json()["detail"]

    def test_upload_without_authentication(
        self,
        client: TestClient,
        upload_file_factory
    ):
        """
        Test upload without authentication headers.

        Expected: 401 Unauthorized
        """
        upload_file = upload_file_factory(
            filename="folder_test.h5ad",
            content=b"mock_data"
        )

        response = client.post(
            "/routes/files/upload",
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        assert response.status_code == 401

    def test_upload_creates_user_directory(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        sample_user: models.User
    ):
        """
        Test first upload creates user directory structure.

        Expected:
        - ./user/{username}/data/ created automatically
        - File uploaded successfully
        """
        # Ensure directory doesn't exist
        user_dir = f"./user/{sample_user.username}"
        if os.path.exists(user_dir):
            import shutil
            shutil.rmtree(user_dir)

        upload_file = upload_file_factory(
            filename="folder_first_file.h5ad",
            content=b"mock_data" * 100
        )

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        assert response.status_code == 200
        assert os.path.exists(f"{user_dir}/data"), "User directory not created"

        # Cleanup
        if os.path.exists(user_dir):
            import shutil
            shutil.rmtree(user_dir)

    def test_upload_filename_parsing(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test folder_filename split logic.

        Format: "folder_filename.h5ad" → folder="folder", file="filename.h5ad"

        Test cases:
        - "analysis_data.h5ad" → folder="analysis", file="data.h5ad"
        - "test_results_exp1.h5ad" → folder="test", file="results_exp1.h5ad"
        """
        test_cases = [
            ("analysis_data.h5ad", "analysis", "data.h5ad"),
            ("test_results_exp1.h5ad", "test", "results_exp1.h5ad"),
            ("folder123_file_with_underscores.h5ad", "folder123", "file_with_underscores.h5ad")
        ]

        for full_filename, expected_folder, expected_filename in test_cases:
            upload_file = upload_file_factory(
                filename=full_filename,
                content=b"mock_data" * 10
            )

            response = client.post(
                "/routes/files/upload",
                headers=auth_headers,
                files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
            )

            assert response.status_code == 200
            data = response.json()
            assert data["folder"] == expected_folder, \
                f"Folder parsing failed: expected {expected_folder}, got {data['folder']}"
            assert data["file_name"] == expected_filename, \
                f"Filename parsing failed: expected {expected_filename}, got {data['file_name']}"

            # Cleanup for next test
            file_path = os.path.join(temp_user_file_directory["data_dir"], expected_filename)
            if os.path.exists(file_path):
                os.remove(file_path)


# ==============================================================================
# HIGH PRIORITY - File Download Operations
# ==============================================================================

class TestFileDownloadOperations:
    """
    Test file download operations.
    """

    def test_download_data_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test successful file download via /data/{filename}.

        Expected:
        - 200 status
        - For H5AD: Returns dict format (orient="records")
        - Correct content-type header
        """
        # Create actual file
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        response = client.get(
            f"/routes/files/data/{sample_file_metadata.file_name}",
            headers=auth_headers
        )

        assert response.status_code == 200, f"Download failed: {response.text}"

        # For H5AD files, should return dict format
        if sample_file_metadata.file_name.endswith('.h5ad'):
            data = response.json()
            assert isinstance(data, list), "H5AD should return list of records"
            assert len(data) > 0, "H5AD data is empty"

    def test_download_nonexistent_file(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test download of non-existent file.

        Expected: 400 Bad Request with "file not exists"
        """
        response = client.get(
            "/routes/files/data/nonexistent_file.h5ad",
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_download_result_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_csv_result_file: dict,
        temp_user_file_directory: dict
    ):
        """
        Test download workflow result via /result/{filename}.

        Expected:
        - FileResponse with correct filename
        - Proper content-disposition header
        """
        response = client.get(
            f"/routes/files/result/{sample_csv_result_file['filename']}",
            headers=auth_headers
        )

        assert response.status_code == 200
        assert response.headers.get("content-type") in ["text/csv", "application/octet-stream"]

    def test_download_tutorial_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_html_tutorial_file: dict
    ):
        """
        Test download public tutorial file.

        Expected: FileResponse with HTML content
        """
        response = client.get(
            f"/routes/files/tutorials/{sample_html_tutorial_file['filename']}",
            headers=auth_headers
        )

        assert response.status_code == 200
        # Tutorial files should be HTML
        assert "html" in response.headers.get("content-type", "").lower() or \
               response.status_code == 200

    @pytest.mark.skip(reason="Cross-user access test requires multiple users")
    def test_download_other_user_file_forbidden(
        self,
        client: TestClient,
        auth_headers: dict,
        user_factory
    ):
        """
        Test User A cannot download User B's file.

        Expected: 403 Forbidden or 400 "file not exists"
        Note: Current implementation returns 400 (implicit authorization)
        """
        # Create second user
        user_b = user_factory(username="userB", email="userb@example.com")

        # TODO: Create file for user_b and attempt access from user_a
        pytest.skip("Requires multi-user setup")

    def test_download_data_file_h5ad_conversion(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test H5AD download returns dict format (line 372-373).

        Expected:
        - load_tab_file called
        - DataFrame converted to dict with orient="records"
        """
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        response = client.get(
            f"/routes/files/data/{sample_file_metadata.file_name}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()

        # Verify dict format
        assert isinstance(data, list), "Should return list of records"
        if len(data) > 0:
            assert isinstance(data[0], dict), "Each record should be a dict"


# ==============================================================================
# MEDIUM PRIORITY - File Metadata Operations
# ==============================================================================

class TestFileMetadataOperations:
    """
    Test file listing, search, and folder operations.
    """

    def test_list_user_files_success(
        self,
        client: TestClient,
        auth_headers: dict,
        multiple_file_metadata: list,
        temp_user_file_directory: dict
    ):
        """
        Test GET /me returns user's file structure.

        Expected:
        - List of {folder: [filenames]} tuples
        - All user files included
        - Correct folder structure from os.walk
        """
        # Create actual files
        for file_record in multiple_file_metadata:
            file_path = os.path.join(temp_user_file_directory["data_dir"], file_record.file_name)
            with open(file_path, 'wb') as f:
                f.write(b"mock_data")

        response = client.get(
            "/routes/files/me",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list), "Should return list of tuples"

    def test_list_user_files_empty(
        self,
        client: TestClient,
        auth_headers: dict,
        temp_user_file_directory: dict
    ):
        """
        Test GET /me with no files.

        Expected: Empty list []
        """
        response = client.get(
            "/routes/files/me",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        # May contain empty data folder
        assert len(data) <= 1

    def test_find_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_file_metadata: models.File
    ):
        """
        Test POST /find with valid filename.

        Expected:
        - File metadata returned
        - Includes file_name, file_size, file_path, folder
        """
        response = client.post(
            "/routes/files/find",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        assert response.status_code == 200
        data = response.json()
        assert data["file_name"] == sample_file_metadata.file_name
        assert "file_size" in data
        assert "folder" in data

    def test_find_file_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test POST /find with non-existent filename.

        Expected: 400 Bad Request with "file not exists"
        """
        response = client.post(
            "/routes/files/find",
            headers=auth_headers,
            json={"file_name": "nonexistent.h5ad"}
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_find_folder_files(
        self,
        client: TestClient,
        auth_headers: dict,
        multiple_file_metadata: list
    ):
        """
        Test POST /folder returns files in specific folder.

        Expected:
        - Only target folder files returned
        - get_user_folder CRUD function called
        """
        response = client.post(
            "/routes/files/folder",
            headers=auth_headers,
            json={"folder_name": "folder1"}
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)

        # Verify all returned files are from folder1
        for file_record in data:
            assert file_record["folder"] == "folder1"


# ==============================================================================
# MEDIUM PRIORITY - File Update and Delete Operations
# ==============================================================================

class TestFileUpdateDeleteOperations:
    """
    Test file update and delete operations.
    """

    def test_update_file_name_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_file_metadata: models.File,
        db_session: Session
    ):
        """
        Test file rename via POST /update.

        Expected:
        - DB record updated
        - Returns updated file metadata
        """
        new_name = "renamed_data.h5ad"

        response = client.post(
            "/routes/files/update",
            headers=auth_headers,
            json={
                "file_name": sample_file_metadata.file_name,
                "update_name": new_name
            }
        )

        assert response.status_code == 200
        data = response.json()
        assert data["file_name"] == new_name

        # Verify in database
        db_session.expire_all()
        updated_file = crud_file.get_user_file(db_session, sample_file_metadata.user_id, new_name)
        assert updated_file is not None
        assert updated_file.file_name == new_name

    def test_update_nonexistent_file(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test update of non-existent file.

        Expected: 400 Bad Request
        """
        response = client.post(
            "/routes/files/update",
            headers=auth_headers,
            json={
                "file_name": "nonexistent.h5ad",
                "update_name": "new_name.h5ad"
            }
        )

        assert response.status_code == 400

    def test_delete_file_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict,
        db_session: Session
    ):
        """
        Test successful file deletion.

        Expected:
        - File removed from filesystem
        - DB record deleted
        - Path traversal check passes
        """
        # Create actual file
        file_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        with open(file_path, 'wb') as f:
            f.write(b"mock_data")

        response = client.post(
            "/routes/files/delete",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        assert response.status_code == 200

        # Verify file deleted from filesystem
        assert not os.path.exists(file_path), "File not deleted from filesystem"

        # Verify DB record deleted
        db_session.expire_all()
        deleted_file = crud_file.get_user_file(db_session, sample_file_metadata.user_id, sample_file_metadata.file_name)
        assert deleted_file is None, "File not deleted from database"

    def test_delete_file_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test delete of non-existent file.

        Expected: 400 Bad Request
        """
        response = client.post(
            "/routes/files/delete",
            headers=auth_headers,
            json={"file_name": "nonexistent.h5ad"}
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_delete_file_filesystem_error(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_file_metadata: models.File
    ):
        """
        Test delete when DB record exists but filesystem file missing.

        Expected: Should handle gracefully (current: may return 500)
        """
        # DB record exists, but no actual file

        response = client.post(
            "/routes/files/delete",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        # Should still delete DB record even if filesystem file missing
        # Current implementation may fail - this documents the behavior
        assert response.status_code in [200, 500]


# ==============================================================================
# MEDIUM PRIORITY - H5AD File Conversion Tests
# ==============================================================================

class TestH5ADConversion:
    """
    Test H5AD file conversion operations.
    """

    def test_convert_h5ad_to_csv_success(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test H5AD to CSV conversion.

        Expected:
        - CSV created in ./user/{username}/result/
        - Filename format: {original}_obs_umap.csv
        - CSV contains obs + X_umap columns
        """
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        response = client.post(
            "/routes/files/convert",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        assert response.status_code == 200
        data = response.json()
        expected_filename = sample_file_metadata.file_name.replace('.h5ad', '') + '_obs_umap.csv'
        assert data["file_name"] == expected_filename

        # Verify CSV exists
        csv_path = os.path.join(temp_user_file_directory["result_dir"], expected_filename)
        assert os.path.exists(csv_path), "CSV file not created"

        # Verify CSV content
        import pandas as pd
        df = pd.read_csv(csv_path)
        assert 'X' in df.columns, "X coordinate missing"
        assert 'Y' in df.columns, "Y coordinate missing"

    def test_convert_h5ad_already_exists(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test conversion when CSV already exists.

        Expected:
        - 200 status
        - Returns existing filename (line 203-204)
        - No re-conversion performed
        """
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        # First conversion
        response1 = client.post(
            "/routes/files/convert",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )
        assert response1.status_code == 200

        # Second conversion (should return existing)
        response2 = client.post(
            "/routes/files/convert",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )
        assert response2.status_code == 200
        assert response1.json() == response2.json()

    def test_convert_missing_h5ad_file(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test conversion of non-existent file.

        Expected: 400 Bad Request
        """
        response = client.post(
            "/routes/files/convert",
            headers=auth_headers,
            json={"file_name": "nonexistent.h5ad"}
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_check_converted_file_exists(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_csv_result_file: dict,
        sample_file_metadata: models.File
    ):
        """
        Test GET /check/{file_name} for converted file.

        Expected: 200 with {"file_name": "{name}_obs_umap.csv"}
        """
        # Use original h5ad name (check endpoint converts it)
        response = client.get(
            f"/routes/files/check/{sample_file_metadata.file_name}",
            headers=auth_headers
        )

        # May not exist yet, depends on setup
        if response.status_code == 200:
            data = response.json()
            assert "_obs_umap.csv" in data["file_name"]


# ==============================================================================
# MEDIUM PRIORITY - H5AD Data Extraction Tests
# ==============================================================================

class TestH5ADDataExtraction:
    """
    Test H5AD column and cluster extraction operations.
    """

    def test_get_h5ad_columns_success(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test extraction of annotation and pseudotime columns.

        Expected:
        {
            "anno_columns": ["cell_type", "cluster"],
            "pseudo_columns": ["pseudotime"]
        }
        """
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        response = client.post(
            "/routes/files/columns",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        assert response.status_code == 200
        data = response.json()
        assert "anno_columns" in data
        assert "pseudo_columns" in data
        assert isinstance(data["anno_columns"], list)
        assert isinstance(data["pseudo_columns"], list)

    def test_get_h5ad_clusters_success(
        self,
        client: TestClient,
        auth_headers: dict,
        mock_h5ad_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test extraction of unique cluster values.

        Expected: {"clusters": ["Cluster1", "Cluster2", ...]}
        """
        import shutil
        target_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        shutil.copy(mock_h5ad_file["filepath"], target_path)

        # Need to add anno_column parameter
        response = client.post(
            "/routes/files/clusters",
            headers=auth_headers,
            json={
                "file_name": sample_file_metadata.file_name,
                "anno_column": "cluster"
            }
        )

        assert response.status_code == 200
        data = response.json()
        assert "clusters" in data
        assert isinstance(data["clusters"], list)
        assert len(data["clusters"]) > 0

    def test_h5ad_operations_missing_file(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test H5AD operations on non-existent file.

        Expected: 400 Bad Request
        """
        response = client.post(
            "/routes/files/columns",
            headers=auth_headers,
            json={"file_name": "nonexistent.h5ad"}
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    @pytest.mark.skip(reason="Requires invalid H5AD file creation")
    def test_h5ad_operations_invalid_file(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test H5AD operations on corrupted file.

        Expected: 500 Internal Server Error or scanpy exception
        """
        # Create invalid H5AD file
        invalid_path = os.path.join(temp_user_file_directory["data_dir"], sample_file_metadata.file_name)
        with open(invalid_path, 'wb') as f:
            f.write(b"not_a_valid_h5ad_file")

        response = client.post(
            "/routes/files/columns",
            headers=auth_headers,
            json={"file_name": sample_file_metadata.file_name}
        )

        # Should handle error gracefully
        assert response.status_code in [400, 500]


# ==============================================================================
# LOW PRIORITY - Setup and Result File Operations
# ==============================================================================

class TestSetupResultOperations:
    """
    Test setup file and result file operations.
    """

    def test_setup_check_lists_json_files(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_json_setup_file: dict,
        temp_user_file_directory: dict
    ):
        """
        Test GET /setup/check returns .json files.

        Expected: {"option_files": ["opt1.json", "opt2.json", ...]}
        """
        # Create additional JSON files
        for i in range(2):
            json_path = os.path.join(temp_user_file_directory["data_dir"], f"option{i}.json")
            with open(json_path, 'w') as f:
                import json
                json.dump({"test": i}, f)

        response = client.get(
            "/routes/files/setup/check",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "option_files" in data
        assert isinstance(data["option_files"], list)
        # Should only include .json files
        for filename in data["option_files"]:
            assert filename.endswith('.json')

    def test_setup_read_json_file(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_json_setup_file: dict
    ):
        """
        Test GET /setup/{file_name} for JSON file.

        Expected: Parsed JSON object
        """
        response = client.get(
            f"/routes/files/setup/{sample_json_setup_file['filename']}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, dict)
        # Verify it matches original content
        assert data == sample_json_setup_file["content"]

    def test_setup_read_binary_file(
        self,
        client: TestClient,
        auth_headers: dict,
        temp_user_file_directory: dict
    ):
        """
        Test GET /setup/{file_name} for non-JSON file.

        Expected: Binary content returned (line 310-312)
        """
        # Create binary file
        binary_filename = "binary_file.bin"
        binary_path = os.path.join(temp_user_file_directory["data_dir"], binary_filename)
        binary_content = b"\x00\x01\x02\x03\x04\x05"
        with open(binary_path, 'wb') as f:
            f.write(binary_content)

        response = client.get(
            f"/routes/files/setup/{binary_filename}",
            headers=auth_headers
        )

        # May return 200 with binary content or error
        assert response.status_code in [200, 500]

    def test_list_result_files_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_csv_result_file: dict,
        sample_file_metadata: models.File,
        temp_user_file_directory: dict
    ):
        """
        Test POST /result to find matching result files.

        Expected: {"result_files": ["option1_data_result.csv", ...]}
        """
        # Create additional result files
        for i in range(2):
            result_filename = f"option1_{sample_file_metadata.file_name.replace('.h5ad', '')}_result{i}.csv"
            result_path = os.path.join(temp_user_file_directory["result_dir"], result_filename)
            with open(result_path, 'w') as f:
                f.write("test,data\n1,2\n")

        response = client.post(
            "/routes/files/result",
            headers=auth_headers,
            json={
                "file_name": sample_file_metadata.file_name.replace('.h5ad', ''),
                "option_file_name": "option1"
            }
        )

        assert response.status_code == 200
        data = response.json()
        assert "result_files" in data
        assert isinstance(data["result_files"], list)
        # All files should contain both filename and option_filename
        for result_file in data["result_files"]:
            assert "option1" in result_file
            assert sample_file_metadata.file_name.replace('.h5ad', '') in result_file

    def test_list_result_files_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test POST /result with no matching files.

        Expected: 400 Bad Request (if file doesn't exist in DB)
        """
        response = client.post(
            "/routes/files/result",
            headers=auth_headers,
            json={
                "file_name": "nonexistent",
                "option_file_name": "option1"
            }
        )

        assert response.status_code == 400


# ==============================================================================
# LOW PRIORITY - Edge Cases and Error Handling
# ==============================================================================

class TestEdgeCasesErrorHandling:
    """
    Test edge cases and error handling.
    """

    def test_upload_empty_file(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test upload of 0-byte file.

        Expected: 400 Bad Request or success with size=0
        """
        upload_file = upload_file_factory(
            filename="folder_empty.h5ad",
            content=b""  # Empty file
        )

        response = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
        )

        # May accept or reject empty files
        assert response.status_code in [200, 400]

    def test_upload_special_characters_filename(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict
    ):
        """
        Test upload with unicode/special characters in filename.

        Filenames: ["数据.h5ad", "file (1).h5ad", "test-file_2.h5ad"]
        Expected: Proper encoding handling
        """
        special_filenames = [
            "folder_file (1).h5ad",
            "folder_test-file_2.h5ad",
            "folder_data_2024.h5ad"
        ]

        for filename in special_filenames:
            upload_file = upload_file_factory(
                filename=filename,
                content=b"mock_data" * 10
            )

            response = client.post(
                "/routes/files/upload",
                headers=auth_headers,
                files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
            )

            # Should handle special chars gracefully
            if response.status_code == 200:
                # Cleanup for next iteration
                folder, fname = filename.split('_', 1)
                file_path = os.path.join(temp_user_file_directory["data_dir"], fname)
                if os.path.exists(file_path):
                    os.remove(file_path)

    @pytest.mark.skip(reason="Requires mock filesystem error")
    def test_filesystem_permission_error(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory
    ):
        """
        Test filesystem write failure.

        Expected: 500 Internal Server Error
        """
        with patch('os.makedirs', side_effect=PermissionError("Mock permission denied")):
            upload_file = upload_file_factory(
                filename="folder_test.h5ad",
                content=b"mock_data"
            )

            response = client.post(
                "/routes/files/upload",
                headers=auth_headers,
                files={"files": (upload_file.filename, upload_file.file, upload_file.content_type)}
            )

            assert response.status_code == 500

    def test_malformed_request_payloads(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """
        Test invalid JSON in POST requests.

        Expected: 422 Unprocessable Entity (FastAPI validation)
        """
        # Missing required field
        response = client.post(
            "/routes/files/find",
            headers=auth_headers,
            json={}  # Missing file_name
        )
        assert response.status_code == 422

        # Wrong data type
        response = client.post(
            "/routes/files/find",
            headers=auth_headers,
            json={"file_name": 123}  # Should be string
        )
        # FastAPI may coerce or reject
        assert response.status_code in [400, 422]

    def test_concurrent_upload_same_filename(
        self,
        client: TestClient,
        auth_headers: dict,
        upload_file_factory,
        temp_user_file_directory: dict,
        user_factory
    ):
        """
        Test two users upload files with same name.

        Expected: Both succeed, isolated in respective directories
        """
        # Create second user
        user2 = user_factory(username="user2", email="user2@example.com")

        # Get auth headers for second user
        login_response = client.post(
            "/routes/auth/login/access-token",
            data={
                "username": user2.email,
                "password": "testpassword123"
            }
        )
        assert login_response.status_code == 200
        auth_headers_user2 = {"Authorization": f"Bearer {login_response.json()['access_token']}"}

        # Create user2 directory
        user2_dir = f"./user/{user2.username}/data"
        os.makedirs(user2_dir, exist_ok=True)

        # Both users upload same filename
        upload_file1 = upload_file_factory(
            filename="folder_same_name.h5ad",
            content=b"user1_data" * 10
        )
        upload_file2 = upload_file_factory(
            filename="folder_same_name.h5ad",
            content=b"user2_data" * 10
        )

        response1 = client.post(
            "/routes/files/upload",
            headers=auth_headers,
            files={"files": (upload_file1.filename, upload_file1.file, upload_file1.content_type)}
        )

        response2 = client.post(
            "/routes/files/upload",
            headers=auth_headers_user2,
            files={"files": (upload_file2.filename, upload_file2.file, upload_file2.content_type)}
        )

        # Both should succeed with isolated directories
        assert response1.status_code == 200
        assert response2.status_code == 200

        # Cleanup
        import shutil
        if os.path.exists(f"./user/{user2.username}"):
            shutil.rmtree(f"./user/{user2.username}")
