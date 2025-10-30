"""
Unit tests for file security utilities.

Tests the validate_file_path() function's ability to prevent
path traversal attacks, including sibling directory attacks.
"""
import os
import pytest
import tempfile
from fastapi import HTTPException
from pathlib import Path

from app.common.utils.file_security import (
    validate_file_path,
    sanitize_filename,
    validate_file_upload
)


class TestValidateFilePath:
    """Test path traversal prevention in validate_file_path()."""

    def test_valid_file_path_allowed(self):
        """Test that valid file paths within allowed directory are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_folder = os.path.join(tmpdir, "user", "data")
            os.makedirs(user_folder, exist_ok=True)

            # Create a valid file
            file_path = os.path.join(user_folder, "test.h5ad")
            Path(file_path).touch()

            # Should not raise exception
            validate_file_path(user_folder, file_path)

    def test_parent_directory_traversal_blocked(self):
        """Test that parent directory traversal (../) is blocked."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_folder = os.path.join(tmpdir, "user", "data")
            os.makedirs(user_folder, exist_ok=True)

            # Attempt to access parent directory
            file_path = os.path.join(user_folder, "../../../etc/passwd")

            with pytest.raises(HTTPException) as exc_info:
                validate_file_path(user_folder, file_path)

            assert exc_info.value.status_code == 403
            assert "Path traversal attempt detected" in exc_info.value.detail

    def test_sibling_directory_attack_blocked(self):
        """
        Test CRITICAL sibling directory attack prevention.

        Attack: If user folder is /tmp/user/john/data
        Attacker provides: ../data_backup/secrets.txt
        After realpath: /tmp/user/john/data_backup/secrets.txt

        Old vulnerable code using startswith():
        - This PASSES because path starts with "/tmp/user/john/data"

        New secure code using Path.relative_to():
        - This FAILS with ValueError because path is not under allowed directory
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Setup: Create user folder and sibling directory
            user_base = os.path.join(tmpdir, "user", "john")
            user_folder = os.path.join(user_base, "data")
            sibling_folder = os.path.join(user_base, "data_backup")

            os.makedirs(user_folder, exist_ok=True)
            os.makedirs(sibling_folder, exist_ok=True)

            # Create "secret" file in sibling directory
            secret_file = os.path.join(sibling_folder, "secrets.txt")
            with open(secret_file, 'w') as f:
                f.write("CONFIDENTIAL DATA")

            # Attempt sibling directory attack
            malicious_path = os.path.join(user_folder, "../data_backup/secrets.txt")

            # Verify file actually exists (to confirm the attack is viable)
            assert os.path.exists(os.path.realpath(malicious_path)), \
                "Test setup failed: secret file should exist"

            # Our security function must block this
            with pytest.raises(HTTPException) as exc_info:
                validate_file_path(user_folder, malicious_path)

            assert exc_info.value.status_code == 403
            assert "Path traversal attempt detected" in exc_info.value.detail

    def test_symlink_attack_blocked(self):
        """Test that symlink attacks are blocked."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_folder = os.path.join(tmpdir, "user", "data")
            os.makedirs(user_folder, exist_ok=True)

            # Create a target file outside user folder
            secret_dir = os.path.join(tmpdir, "secrets")
            os.makedirs(secret_dir, exist_ok=True)
            secret_file = os.path.join(secret_dir, "password.txt")
            with open(secret_file, 'w') as f:
                f.write("SECRET")

            # Create symlink inside user folder pointing to secret file
            symlink_path = os.path.join(user_folder, "innocent.txt")
            try:
                os.symlink(secret_file, symlink_path)
            except OSError:
                pytest.skip("Symlink creation not supported on this system")

            # Attempt to access via symlink
            with pytest.raises(HTTPException) as exc_info:
                validate_file_path(user_folder, symlink_path)

            assert exc_info.value.status_code == 403
            assert "Path traversal attempt detected" in exc_info.value.detail

    def test_absolute_path_outside_folder_blocked(self):
        """Test that absolute paths outside allowed folder are blocked."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_folder = os.path.join(tmpdir, "user", "data")
            os.makedirs(user_folder, exist_ok=True)

            # Attempt to access absolute path outside user folder
            malicious_path = "/etc/passwd"

            with pytest.raises(HTTPException) as exc_info:
                validate_file_path(user_folder, malicious_path)

            assert exc_info.value.status_code == 403
            assert "Path traversal attempt detected" in exc_info.value.detail

    def test_nested_subdirectory_allowed(self):
        """Test that nested subdirectories within user folder are allowed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_folder = os.path.join(tmpdir, "user", "data")
            nested_dir = os.path.join(user_folder, "subdir", "nested")
            os.makedirs(nested_dir, exist_ok=True)

            file_path = os.path.join(nested_dir, "test.h5ad")
            Path(file_path).touch()

            # Should not raise exception
            validate_file_path(user_folder, file_path)

    def test_same_prefix_different_folder_blocked(self):
        """Test that folders with same prefix but different paths are blocked."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Setup: Create two folders with similar names
            user_folder = os.path.join(tmpdir, "user", "data")
            other_folder = os.path.join(tmpdir, "user", "data_other")

            os.makedirs(user_folder, exist_ok=True)
            os.makedirs(other_folder, exist_ok=True)

            # Create file in other folder
            other_file = os.path.join(other_folder, "file.txt")
            Path(other_file).touch()

            # Attempt to access via path traversal
            malicious_path = os.path.join(user_folder, "../data_other/file.txt")

            with pytest.raises(HTTPException) as exc_info:
                validate_file_path(user_folder, malicious_path)

            assert exc_info.value.status_code == 403
            assert "Path traversal attempt detected" in exc_info.value.detail


class TestSanitizeFilename:
    """Test filename sanitization."""

    def test_sanitize_dangerous_characters(self):
        """Test that dangerous characters are removed."""
        # Spaces and special chars are replaced with underscores
        assert sanitize_filename("file;rm -rf /.h5ad") == "file_rm_-rf__.h5ad"
        assert sanitize_filename("file<script>.h5ad") == "file_script_.h5ad"
        assert sanitize_filename("file\x00.h5ad") == "file_.h5ad"

    def test_sanitize_path_traversal_characters(self):
        """Test that path traversal characters are sanitized."""
        # Slashes are replaced with underscores
        assert sanitize_filename("../../etc/passwd") == ".._.._etc_passwd"
        assert sanitize_filename("../secrets.txt") == ".._secrets.txt"

    def test_sanitize_long_filename(self):
        """Test that long filenames are truncated."""
        long_name = "a" * 300 + ".h5ad"
        sanitized = sanitize_filename(long_name, max_length=255)
        assert len(sanitized) <= 255
        assert sanitized.endswith(".h5ad")

    def test_sanitize_preserves_valid_filenames(self):
        """Test that valid filenames are preserved."""
        assert sanitize_filename("valid_file-123.h5ad") == "valid_file-123.h5ad"
        assert sanitize_filename("file.with.dots.csv") == "file.with.dots.csv"


class TestValidateFileUploadWithMIME:
    """Test enhanced file upload validation with MIME type checking."""

    class MockUploadFile:
        """Mock FastAPI UploadFile for testing."""

        def __init__(self, filename: str, content: bytes, content_type: str = "application/octet-stream"):
            self.filename = filename
            self.content_type = content_type
            from io import BytesIO
            self.file = BytesIO(content)

    def test_validate_with_mime_check_valid_h5ad(self, monkeypatch):
        """Test validation with valid H5AD file and MIME checking."""
        monkeypatch.setattr("app.common.config.ENABLE_FILE_SIGNATURE_VALIDATION", True)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create valid H5AD file
            temp_path = os.path.join(tmpdir, "temp.h5ad")
            with open(temp_path, 'wb') as f:
                f.write(b'\x89HDF\r\n\x1a\n')  # HDF5 signature
                f.write(b'\x00' * 1000)

            # Create mock upload file
            mock_file = self.MockUploadFile(
                filename="test.h5ad",
                content=b'\x89HDF\r\n\x1a\n' + b'\x00' * 1000
            )

            # Should pass validation
            safe_name = validate_file_upload(
                mock_file,
                user_id=123,
                temp_path=temp_path
            )
            assert safe_name == "test.h5ad"

    def test_validate_with_mime_check_spoofed_file(self, monkeypatch, temp_security_log):
        """Test that spoofed files are rejected with MIME checking."""
        monkeypatch.setattr("app.common.config.ENABLE_FILE_SIGNATURE_VALIDATION", True)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create file with wrong signature (executable as H5AD)
            temp_path = os.path.join(tmpdir, "temp.h5ad")
            with open(temp_path, 'wb') as f:
                f.write(b'MZ\x90\x00')  # EXE signature
                f.write(b'\x00' * 1000)

            # Create mock upload file
            mock_file = self.MockUploadFile(
                filename="malware.h5ad",
                content=b'MZ\x90\x00' + b'\x00' * 1000
            )

            # Should raise exception due to signature mismatch
            with pytest.raises(HTTPException) as exc_info:
                validate_file_upload(
                    mock_file,
                    user_id=123,
                    temp_path=temp_path
                )

            assert exc_info.value.status_code == 400
            assert "does not match declared extension" in exc_info.value.detail

    def test_validate_logs_security_events(self, monkeypatch, temp_security_log):
        """Test that security events are logged during validation."""
        # Create mock upload file with invalid extension
        mock_file = self.MockUploadFile(
            filename="test.exe",
            content=b'\x00' * 1000
        )

        # Should raise exception and log event
        with pytest.raises(HTTPException):
            validate_file_upload(mock_file, user_id=456)

        # Verify security log was created
        from app.common.utils.security_logging import get_security_events
        events = get_security_events(event_type="invalid_file_extension")
        assert len(events) > 0
        assert events[0]["user_id"] == 456

    def test_validate_logs_file_size_exceeded(self, monkeypatch, temp_security_log):
        """Test that file size violations are logged."""
        # Create oversized file
        large_content = b'\x00' * (600 * 1024 * 1024)  # 600MB
        mock_file = self.MockUploadFile(
            filename="large.h5ad",
            content=large_content
        )

        # Should raise exception and log event
        with pytest.raises(HTTPException) as exc_info:
            validate_file_upload(mock_file, user_id=789)

        assert exc_info.value.status_code == 413

        # Verify security log
        from app.common.utils.security_logging import get_security_events
        events = get_security_events(event_type="file_size_exceeded")
        assert len(events) > 0
        assert events[0]["user_id"] == 789

    def test_validate_logs_empty_file(self, monkeypatch, temp_security_log):
        """Test that empty file uploads are logged."""
        mock_file = self.MockUploadFile(
            filename="empty.h5ad",
            content=b''
        )

        # Should raise exception and log event
        with pytest.raises(HTTPException):
            validate_file_upload(mock_file, user_id=111)

        # Verify security log
        from app.common.utils.security_logging import get_security_events
        events = get_security_events(event_type="empty_file_upload")
        assert len(events) > 0

    def test_mime_validation_disabled(self, monkeypatch):
        """Test that MIME validation can be disabled via config."""
        monkeypatch.setattr("app.common.config.ENABLE_FILE_SIGNATURE_VALIDATION", False)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create file with wrong signature
            temp_path = os.path.join(tmpdir, "temp.h5ad")
            with open(temp_path, 'wb') as f:
                f.write(b'MZ\x90\x00')  # Wrong signature
                f.write(b'\x00' * 1000)

            mock_file = self.MockUploadFile(
                filename="test.h5ad",
                content=b'MZ\x90\x00' + b'\x00' * 1000
            )

            # Should pass when MIME validation is disabled
            safe_name = validate_file_upload(
                mock_file,
                user_id=123,
                temp_path=temp_path
            )
            assert safe_name == "test.h5ad"


@pytest.fixture
def temp_security_log(monkeypatch):
    """Create a temporary security log file for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        log_file = os.path.join(tmpdir, "test_security.log")
        monkeypatch.setattr("app.common.config.SECURITY_LOG_FILE", log_file)

        # Reset logger handlers
        import logging
        logger = logging.getLogger("security")
        logger.handlers.clear()

        yield log_file

        logger.handlers.clear()
