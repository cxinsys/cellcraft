"""
File security utilities for path validation and sanitization.

This module provides security utilities to prevent path traversal attacks
and filename injection vulnerabilities in file operations.
"""
import os
import re
from pathlib import Path
from typing import Optional
from fastapi import HTTPException


# File signature mappings (magic numbers for content-based validation)
# These are the first few bytes that identify file types reliably
FILE_SIGNATURES = {
    '.h5ad': [
        b'\x89HDF\r\n\x1a\n',  # HDF5 format (H5AD is built on HDF5)
    ],
    '.csv': [
        # CSV has no magic number, validated by structure
    ],
    '.json': [
        # JSON files typically start with { or [ (whitespace allowed)
    ],
    '.txt': [
        # Plain text has no specific magic number
    ],
}


def validate_file_signature(file_path: str, expected_extension: str) -> bool:
    """
    Validate file content matches expected type using magic numbers.

    Prevents extension spoofing where malicious files (e.g., executables)
    are renamed to appear as safe data files. This is a critical defense
    against OWASP A04:2021 - Insecure Design attacks.

    Args:
        file_path: Path to file for validation
        expected_extension: Expected file extension (e.g., '.h5ad', '.csv')

    Returns:
        True if file signature matches expected type, False otherwise

    Raises:
        HTTPException(400): If file cannot be read or signature check fails

    Example:
        >>> validate_file_signature('/tmp/data.h5ad', '.h5ad')  # True for real H5AD
        >>> validate_file_signature('/tmp/malware.exe.h5ad', '.h5ad')  # False for spoofed

    Note:
        This function reads the first 2048 bytes to detect file signatures
        without loading entire files into memory. For formats without magic
        numbers (CSV, JSON, TXT), we perform basic structure validation.
    """
    try:
        # Read first 2048 bytes for signature detection
        with open(file_path, 'rb') as f:
            file_header = f.read(2048)

        if not file_header:
            return False

        # Get expected signatures for this file type
        signatures = FILE_SIGNATURES.get(expected_extension.lower(), [])

        # HDF5/H5AD files have a specific magic number
        if expected_extension.lower() == '.h5ad':
            # Check for HDF5 signature at the beginning
            if not file_header.startswith(b'\x89HDF\r\n\x1a\n'):
                return False
            return True

        # CSV validation: Check for text content with typical CSV structure
        elif expected_extension.lower() == '.csv':
            try:
                # Decode as text and check for CSV-like structure
                text_content = file_header.decode('utf-8', errors='strict')

                # Basic heuristics: should have printable characters and commas
                if ',' in text_content or '\t' in text_content:
                    # Check that content is mostly printable
                    printable_ratio = sum(c.isprintable() or c in '\n\r\t'
                                        for c in text_content) / len(text_content)
                    return printable_ratio > 0.95
                return False
            except UnicodeDecodeError:
                # Binary content is not a valid CSV
                return False

        # JSON validation: Check for valid JSON structure
        elif expected_extension.lower() == '.json':
            try:
                # Decode and check for JSON structure
                text_content = file_header.decode('utf-8', errors='strict').lstrip()

                # JSON must start with { or [
                if text_content and text_content[0] in ('{', '['):
                    # Additional validation: check for mostly printable characters
                    printable_ratio = sum(c.isprintable() or c in '\n\r\t'
                                        for c in text_content) / len(text_content)
                    return printable_ratio > 0.95
                return False
            except UnicodeDecodeError:
                # Binary content is not valid JSON
                return False

        # TXT validation: Check for text content
        elif expected_extension.lower() == '.txt':
            try:
                # Decode as text
                text_content = file_header.decode('utf-8', errors='strict')

                # Check that content is mostly printable
                printable_ratio = sum(c.isprintable() or c in '\n\r\t\x0c'
                                    for c in text_content) / len(text_content)
                return printable_ratio > 0.90  # Allow slightly more binary for TXT
            except UnicodeDecodeError:
                # Try common encodings before rejecting
                for encoding in ['latin-1', 'cp1252']:
                    try:
                        text_content = file_header.decode(encoding)
                        printable_ratio = sum(c.isprintable() or c in '\n\r\t\x0c'
                                            for c in text_content) / len(text_content)
                        if printable_ratio > 0.90:
                            return True
                    except:
                        continue
                return False

        # Unknown extension - allow by default (handled by extension whitelist)
        return True

    except IOError:
        # File read error
        return False


def validate_file_path(user_folder: str, file_path: str, user_id: Optional[int] = None) -> None:
    """
    Validate file path to prevent path traversal attacks.

    Uses pathlib.Path.relative_to() to ensure the resolved file path
    is strictly within the allowed user folder. This prevents sibling
    directory attacks that can bypass startswith() checks.

    Args:
        user_folder: Allowed base directory (e.g., './user/username/data')
        file_path: Full path to file being accessed
        user_id: Optional user ID for security logging

    Raises:
        HTTPException(403): If path traversal detected

    Example:
        validate_file_path('./user/john/data', './user/john/data/file.h5ad')  # OK
        validate_file_path('./user/john/data', './user/john/data/../../../etc/passwd')  # Raises 403
        validate_file_path('./user/john/data', '../data_backup/secret.txt')  # Raises 403 (sibling attack)
    """
    try:
        # Resolve to absolute paths (handles symlinks and relative paths)
        resolved_folder = Path(user_folder).resolve()
        resolved_file = Path(file_path).resolve()

        # Verify file is within allowed folder using relative_to()
        # This raises ValueError if resolved_file is not under resolved_folder
        resolved_file.relative_to(resolved_folder)
    except ValueError:
        # Path is outside allowed directory - log security event
        from app.core.logging import log_security_event

        log_security_event(
            event_type="path_traversal_attempt",
            user_id=user_id,
            details={
                "attempted_path": file_path,
                "user_folder": user_folder,
                "resolved_path": str(Path(file_path).resolve()) if Path(file_path).exists() else None
            },
            severity="CRITICAL"
        )

        raise HTTPException(
            status_code=403,
            detail="Access denied: Path traversal attempt detected"
        )


def sanitize_filename(filename: str, max_length: int = 255) -> str:
    """
    Sanitize filename to prevent injection attacks.

    - Removes special characters (keeps alphanumeric, dash, underscore, dot)
    - Limits filename length
    - Preserves file extension

    Args:
        filename: Original filename
        max_length: Maximum filename length

    Returns:
        Sanitized filename

    Example:
        sanitize_filename("file;rm -rf /.h5ad") -> "filerm_-rf_.h5ad"
        sanitize_filename("../../etc/passwd") -> "..etc_passwd"
    """
    # Replace dangerous characters with underscore
    safe_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', filename)

    # Truncate if too long (preserve extension)
    if len(safe_name) > max_length:
        name, ext = os.path.splitext(safe_name)
        name = name[:max_length - len(ext)]
        safe_name = name + ext

    return safe_name


def validate_file_upload(
    file,
    max_size: int = None,
    allowed_extensions: set = None,
    max_filename_length: int = None,
    user_id: Optional[int] = None,
    temp_path: Optional[str] = None
) -> str:
    """
    Comprehensive file upload validation.

    Validates:
    1. Filename length
    2. File extension (whitelist)
    3. File size
    4. Filename sanitization
    5. MIME type / file signature (if temp_path provided)

    Args:
        file: FastAPI UploadFile object
        max_size: Maximum file size in bytes (default from config)
        allowed_extensions: Set of allowed file extensions (default from config)
        max_filename_length: Maximum filename length (default from config)
        user_id: Optional user ID for security logging
        temp_path: Optional path to temporary file for MIME validation

    Returns:
        Sanitized filename

    Raises:
        HTTPException(400): Invalid filename, extension, size, or MIME type
        HTTPException(413): File too large

    Example:
        safe_name = validate_file_upload(upload_file, user_id=123, temp_path="/tmp/upload")
        # Raises exception if invalid, returns safe filename if valid
    """
    from app.core.config import (
        MAX_UPLOAD_SIZE,
        ALLOWED_EXTENSIONS,
        MAX_FILENAME_LENGTH,
        ENABLE_FILE_SIGNATURE_VALIDATION
    )
    from app.core.logging import log_security_event

    # Use defaults from config if not provided
    if max_size is None:
        max_size = MAX_UPLOAD_SIZE
    if allowed_extensions is None:
        allowed_extensions = ALLOWED_EXTENSIONS
    if max_filename_length is None:
        max_filename_length = MAX_FILENAME_LENGTH

    filename = file.filename

    # 1. Check filename length
    if len(filename) > max_filename_length:
        log_security_event(
            event_type="filename_too_long",
            user_id=user_id,
            details={
                "filename": filename,
                "length": len(filename),
                "max_length": max_filename_length
            },
            severity="WARNING"
        )
        raise HTTPException(
            status_code=400,
            detail=f"Filename too long. Maximum {max_filename_length} characters."
        )

    # 2. Check file extension (whitelist)
    ext = os.path.splitext(filename)[1].lower()
    if ext not in allowed_extensions:
        log_security_event(
            event_type="invalid_file_extension",
            user_id=user_id,
            details={
                "filename": filename,
                "extension": ext,
                "allowed_extensions": list(allowed_extensions)
            },
            severity="WARNING"
        )
        raise HTTPException(
            status_code=400,
            detail=f"Invalid file type '{ext}'. Allowed: {', '.join(sorted(allowed_extensions))}"
        )

    # 3. Check file size (use seek to get size without reading entire file)
    file.file.seek(0, 2)  # Seek to end
    size = file.file.tell()
    file.file.seek(0)  # Reset to beginning

    if size > max_size:
        max_size_mb = max_size / (1024 * 1024)
        actual_size_mb = size / (1024 * 1024)
        log_security_event(
            event_type="file_size_exceeded",
            user_id=user_id,
            details={
                "filename": filename,
                "size_bytes": size,
                "size_mb": actual_size_mb,
                "max_size_mb": max_size_mb
            },
            severity="WARNING"
        )
        raise HTTPException(
            status_code=413,
            detail=f"File too large ({actual_size_mb:.1f}MB). Maximum allowed: {max_size_mb:.0f}MB"
        )

    if size == 0:
        log_security_event(
            event_type="empty_file_upload",
            user_id=user_id,
            details={"filename": filename},
            severity="WARNING"
        )
        raise HTTPException(
            status_code=400,
            detail="Empty file not allowed"
        )

    # 4. Sanitize filename
    safe_filename = sanitize_filename(filename, max_filename_length)

    # 5. MIME type validation (if enabled and temp file provided)
    if ENABLE_FILE_SIGNATURE_VALIDATION and temp_path and os.path.exists(temp_path):
        if not validate_file_signature(temp_path, ext):
            log_security_event(
                event_type="file_signature_mismatch",
                user_id=user_id,
                details={
                    "filename": filename,
                    "declared_extension": ext,
                    "temp_path": temp_path
                },
                severity="CRITICAL"
            )
            raise HTTPException(
                status_code=400,
                detail=f"File content does not match declared extension '{ext}'. Possible file spoofing attempt."
            )

    return safe_filename
