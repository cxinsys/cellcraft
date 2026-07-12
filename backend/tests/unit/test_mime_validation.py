"""
Unit tests for MIME type validation and file signature checking.

Tests the validate_file_signature() function's ability to detect
file content spoofing and prevent malicious file uploads.
"""
import os
import pytest
import tempfile
from pathlib import Path

from app.file.security import validate_file_signature


class TestValidateFileSignature:
    """Test MIME type validation and file signature checking."""

    def test_valid_h5ad_signature(self):
        """Test that valid H5AD files are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a file with valid HDF5 signature
            file_path = os.path.join(tmpdir, "test.h5ad")
            with open(file_path, 'wb') as f:
                # Write HDF5 magic number
                f.write(b'\x89HDF\r\n\x1a\n')
                # Add some dummy data
                f.write(b'\x00' * 1000)

            # Should validate successfully
            assert validate_file_signature(file_path, '.h5ad') is True

    def test_invalid_h5ad_signature(self):
        """Test that files with wrong signature are rejected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a file with wrong signature (e.g., executable)
            file_path = os.path.join(tmpdir, "malware.h5ad")
            with open(file_path, 'wb') as f:
                # Write EXE magic number (MZ header)
                f.write(b'MZ\x90\x00')
                f.write(b'\x00' * 1000)

            # Should fail validation
            assert validate_file_signature(file_path, '.h5ad') is False

    def test_empty_file_rejected(self):
        """Test that empty files are rejected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "empty.h5ad")
            Path(file_path).touch()

            # Should fail validation
            assert validate_file_signature(file_path, '.h5ad') is False

    def test_valid_csv_file(self):
        """Test that valid CSV files are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.csv")
            with open(file_path, 'w') as f:
                f.write("name,age,city\n")
                f.write("Alice,30,NYC\n")
                f.write("Bob,25,LA\n")

            # Should validate successfully
            assert validate_file_signature(file_path, '.csv') is True

    def test_csv_with_tabs(self):
        """Test that tab-separated files are accepted as CSV."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.csv")
            with open(file_path, 'w') as f:
                f.write("name\tage\tcity\n")
                f.write("Alice\t30\tNYC\n")

            # Should validate successfully
            assert validate_file_signature(file_path, '.csv') is True

    def test_binary_file_as_csv_rejected(self):
        """Test that binary files posing as CSV are rejected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "malware.csv")
            with open(file_path, 'wb') as f:
                # Write binary data
                f.write(b'\x89PNG\r\n\x1a\n')
                f.write(b'\x00\x01\x02\x03' * 100)

            # Should fail validation
            assert validate_file_signature(file_path, '.csv') is False

    def test_valid_json_file(self):
        """Test that valid JSON files are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.json")
            with open(file_path, 'w') as f:
                f.write('{"name": "test", "value": 123}')

            # Should validate successfully
            assert validate_file_signature(file_path, '.json') is True

    def test_json_array(self):
        """Test that JSON arrays are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.json")
            with open(file_path, 'w') as f:
                f.write('[1, 2, 3, 4, 5]')

            # Should validate successfully
            assert validate_file_signature(file_path, '.json') is True

    def test_json_with_whitespace(self):
        """Test that JSON with leading whitespace is accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.json")
            with open(file_path, 'w') as f:
                f.write('  \n  {"name": "test"}')

            # Should validate successfully
            assert validate_file_signature(file_path, '.json') is True

    def test_invalid_json_structure(self):
        """Test that files not starting with { or [ are rejected as JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "invalid.json")
            with open(file_path, 'w') as f:
                f.write('This is not JSON')

            # Should fail validation
            assert validate_file_signature(file_path, '.json') is False

    def test_binary_as_json_rejected(self):
        """Test that binary files posing as JSON are rejected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "malware.json")
            with open(file_path, 'wb') as f:
                f.write(b'\x00\x01\x02\x03' * 100)

            # Should fail validation
            assert validate_file_signature(file_path, '.json') is False

    def test_valid_txt_file(self):
        """Test that valid text files are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.txt")
            with open(file_path, 'w') as f:
                f.write("This is a test text file.\n")
                f.write("It contains multiple lines.\n")

            # Should validate successfully
            assert validate_file_signature(file_path, '.txt') is True

    def test_txt_with_special_characters(self):
        """Test that text files with special characters are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.txt")
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write("Unicode: αβγδ €£¥\n")
                f.write("Special: \t\n\r\n")

            # Should validate successfully
            assert validate_file_signature(file_path, '.txt') is True

    def test_mostly_binary_as_txt_rejected(self):
        """Test that files with mostly binary content are rejected as TXT."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "binary.txt")
            with open(file_path, 'wb') as f:
                # Write mostly binary data (>10% non-printable)
                f.write(b'\x00\x01\x02\x03' * 100)
                f.write(b'some text')

            # Should fail validation
            assert validate_file_signature(file_path, '.txt') is False

    def test_latin1_encoded_txt(self):
        """Test that Latin-1 encoded text files are accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "latin1.txt")
            with open(file_path, 'w', encoding='latin-1') as f:
                f.write("Café résumé naïve")

            # Should validate successfully
            assert validate_file_signature(file_path, '.txt') is True

    def test_nonexistent_file(self):
        """Test that nonexistent files return False."""
        result = validate_file_signature('/nonexistent/path/file.h5ad', '.h5ad')
        assert result is False

    def test_unknown_extension_allowed(self):
        """Test that unknown extensions are allowed by default."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.unknown")
            with open(file_path, 'w') as f:
                f.write("test content")

            # Should allow unknown extensions (handled by whitelist)
            assert validate_file_signature(file_path, '.unknown') is True

    def test_case_insensitive_extensions(self):
        """Test that extension matching is case-insensitive."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.H5AD")
            with open(file_path, 'wb') as f:
                f.write(b'\x89HDF\r\n\x1a\n')
                f.write(b'\x00' * 100)

            # Should work with uppercase extension
            assert validate_file_signature(file_path, '.H5AD') is True
            assert validate_file_signature(file_path, '.h5ad') is True

    def test_executable_spoofing_detected(self):
        """
        Test CRITICAL: Detect executable files renamed to safe extensions.

        This is the primary security use case for MIME validation.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Windows PE executable magic number
            exe_path = os.path.join(tmpdir, "virus.h5ad")
            with open(exe_path, 'wb') as f:
                f.write(b'MZ')  # DOS/PE executable signature
                f.write(b'\x90\x00\x03\x00\x00\x00')
                f.write(b'\x00' * 1000)

            # Should detect and reject
            assert validate_file_signature(exe_path, '.h5ad') is False

            # ELF executable magic number
            elf_path = os.path.join(tmpdir, "malware.csv")
            with open(elf_path, 'wb') as f:
                f.write(b'\x7fELF')  # ELF executable signature
                f.write(b'\x00' * 1000)

            # Should detect and reject
            assert validate_file_signature(elf_path, '.csv') is False
