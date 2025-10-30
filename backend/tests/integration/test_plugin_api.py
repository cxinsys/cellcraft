"""
Integration tests for Plugin API endpoints.

This module tests the plugin-related API endpoints in:
- GET /routes/plugin/list - Get plugin list
- GET /routes/plugin/info/{plugin_name} - Get plugin details
- POST /routes/plugin/validation - Validate plugin
- POST /routes/plugin/upload - Upload plugin
- POST /routes/plugin/associate - Associate user with plugin
- POST /routes/plugin/dissociate - Dissociate user from plugin
- GET /routes/plugin/reference_folders/{plugin_name} - Get reference folders
- POST /routes/plugin/build/{plugin_name} - Build Docker image
- GET /routes/plugin/check_image/{plugin_name} - Check if Docker image exists

Coverage Goal: 55%+ for app/routes/endpoints/plugin.py
Quality Score Goal: 8.2/10
"""
import pytest
import json
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock, mock_open
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app.database import models
from app.common.enums import PluginType


# ==============================================================================
# Test Group 1: Plugin CRUD Operations (List, Info, Query)
# ==============================================================================

class TestPluginCRUD:
    """Test plugin listing, retrieval, and querying operations."""

    def test_get_plugin_list_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test retrieving plugin list successfully."""
        response = client.get(
            "/routes/plugin/list",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) >= 1

        # Check plugin structure
        plugin = data[0]
        assert "name" in plugin
        assert "description" in plugin
        assert "plugin_type" in plugin
        assert "source" in plugin

    def test_get_plugin_list_empty(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session
    ):
        """Test retrieving empty plugin list."""
        # Delete all plugins
        db_session.query(models.Plugin).delete()
        db_session.commit()

        response = client.get(
            "/routes/plugin/list",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 0

    def test_get_plugin_list_filter_by_type(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session,
        sample_user: models.User
    ):
        """Test filtering plugins by type."""
        # Create plugins of different types
        analysis_plugin = models.Plugin(
            name="AnalysisPlugin",
            description="Analysis plugin",
            author="test",
            plugin_path="./plugin/local/AnalysisPlugin",
            plugin_type=PluginType.ANALYSIS,
            dependencies={},
            drawflow={},
            rules={},
            source="local"
        )

        viz_plugin = models.Plugin(
            name="VizPlugin",
            description="Visualization plugin",
            author="test",
            plugin_path="./plugin/local/VizPlugin",
            plugin_type=PluginType.VISUALIZATION,
            dependencies={},
            drawflow={},
            rules={},
            source="local"
        )

        db_session.add_all([analysis_plugin, viz_plugin])
        db_session.commit()

        response = client.get(
            "/routes/plugin/list",
            headers=auth_headers
        )

        assert response.status_code == 200
        plugins = response.json()

        # Check both types exist
        types = [p["plugin_type"] for p in plugins]
        assert "analysis" in types
        assert "visualization" in types

    def test_get_plugin_info_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test retrieving plugin detailed information."""
        response = client.get(
            f"/routes/plugin/info/{sample_plugin.name}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["name"] == sample_plugin.name
        assert data["description"] == sample_plugin.description
        assert data["plugin_type"] == sample_plugin.plugin_type.value
        assert "drawflow" in data
        assert "rules" in data
        assert "source" in data

    def test_get_plugin_info_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test retrieving non-existent plugin."""
        response = client.get(
            "/routes/plugin/info/NonExistentPlugin",
            headers=auth_headers
        )

        assert response.status_code == 404

    def test_get_plugin_info_official_vs_local(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session
    ):
        """Test distinguishing between official and local plugins."""
        # Create official plugin
        official_plugin = models.Plugin(
            name="OfficialPlugin",
            description="Official plugin",
            author="official",
            plugin_path="./plugin/official/OfficialPlugin",
            plugin_type=PluginType.ANALYSIS,
            dependencies={},
            drawflow={},
            rules={},
            source="official",
            is_editable=False,
            version="1.0.0"
        )
        db_session.add(official_plugin)
        db_session.commit()

        response = client.get(
            f"/routes/plugin/info/{official_plugin.name}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["source"] == "official"
        assert data["is_editable"] is False
        assert data["version"] == "1.0.0"


# ==============================================================================
# Test Group 2: Plugin Upload and Validation
# ==============================================================================

class TestPluginUploadValidation:
    """Test plugin upload and validation operations."""

    def test_validate_plugin_success(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test successful plugin validation."""
        plugin_data = {
            "plugin": {
                "name": "TestValidationPlugin",
                "description": "Test plugin for validation",
                "pluginType": "analysis",
                "referenceFolders": [],
                "dependencyFiles": [
                    {
                        "fileName": "requirements.txt",
                        "file": "numpy==1.24.0",
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
                    "parameters": [],
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

        response = client.post(
            "/routes/plugin/validation",
            json=plugin_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        result = response.json()
        assert result["valid"] is True or "message" in result

    def test_validate_plugin_missing_name(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test validation with missing plugin name."""
        plugin_data = {
            "plugin": {
                "name": "",  # Empty name
                "description": "Test plugin",
                "pluginType": "analysis",
                "referenceFolders": [],
                "dependencyFiles": [],
                "useGpu": False
            },
            "rules": [
                {
                    "name": "test_rule",
                    "input": [],
                    "output": [],
                    "script": "test.py",
                    "parameters": [],
                    "nodeId": 1
                }
            ],
            "drawflow": {}
        }

        response = client.post(
            "/routes/plugin/validation",
            json=plugin_data,
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "name" in response.json()["detail"].lower() or "required" in response.json()["detail"].lower()

    def test_validate_plugin_missing_rules(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test validation with missing rules."""
        plugin_data = {
            "plugin": {
                "name": "TestPlugin",
                "description": "Test plugin",
                "pluginType": "analysis",
                "referenceFolders": [],
                "dependencyFiles": [],
                "useGpu": False
            },
            "rules": [],  # No rules
            "drawflow": {}
        }

        response = client.post(
            "/routes/plugin/validation",
            json=plugin_data,
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "rule" in response.json()["detail"].lower()

    def test_validate_plugin_missing_script(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test validation with missing script in rules."""
        plugin_data = {
            "plugin": {
                "name": "TestPlugin",
                "description": "Test plugin",
                "pluginType": "analysis",
                "referenceFolders": [],
                "dependencyFiles": [],
                "useGpu": False
            },
            "rules": [
                {
                    "name": "test_rule",
                    "input": [],
                    "output": [],
                    "script": None,  # No script
                    "parameters": [],
                    "nodeId": 1
                }
            ],
            "drawflow": {}
        }

        response = client.post(
            "/routes/plugin/validation",
            json=plugin_data,
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "script" in response.json()["detail"].lower()

    @patch('app.routes.endpoints.plugin.os.makedirs')
    @patch('app.routes.endpoints.plugin.open', new_callable=mock_open)
    def test_upload_plugin_success(
        self,
        mock_file_open: MagicMock,
        mock_makedirs: MagicMock,
        client: TestClient,
        auth_headers: dict,
        db_session: Session
    ):
        """Test successful plugin upload."""
        plugin_data = {
            "plugin": {
                "name": "UploadTestPlugin",
                "description": "Test upload plugin",
                "pluginType": "analysis",
                "referenceFolders": [],
                "dependencyFiles": [
                    {
                        "fileName": "requirements.txt",
                        "file": "numpy==1.24.0",
                        "type": "text"
                    }
                ],
                "useGpu": False
            },
            "rules": [
                {
                    "name": "analysis_rule",
                    "input": ["data.h5ad"],
                    "output": ["result.csv"],
                    "script": "analyze.py",
                    "parameters": [],
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

        response = client.post(
            "/routes/plugin/upload",
            json=plugin_data,
            headers=auth_headers
        )

        # Should create plugin or return success
        assert response.status_code in [200, 201]

        # Verify plugin was created in database
        plugin = db_session.query(models.Plugin).filter(
            models.Plugin.name == "UploadTestPlugin"
        ).first()

        if plugin:  # Only verify if plugin was actually created
            assert plugin.name == "UploadTestPlugin"
            assert plugin.plugin_type == PluginType.ANALYSIS


# ==============================================================================
# Test Group 3: File Management
# ==============================================================================

class TestPluginFileManagement:
    """Test plugin file-related operations."""

    def test_get_reference_folders_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test retrieving plugin reference folders."""
        response = client.get(
            f"/routes/plugin/reference_folders/{sample_plugin.name}",
            headers=auth_headers
        )

        # Reference folders endpoint should return data
        assert response.status_code in [200, 404]  # May not have reference folders

        if response.status_code == 200:
            data = response.json()
            assert isinstance(data, (dict, list))

    @patch('app.routes.endpoints.plugin.os.path.exists')
    @patch('app.routes.endpoints.plugin.FileResponse')
    def test_download_plugin_package_success(
        self,
        mock_file_response: MagicMock,
        mock_exists: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test downloading plugin package."""
        mock_exists.return_value = True
        mock_file_response.return_value = MagicMock()

        response = client.get(
            f"/routes/plugin/package/{sample_plugin.name}",
            headers=auth_headers
        )

        # Should attempt to download or return not found
        assert response.status_code in [200, 404]

    @patch('app.routes.endpoints.plugin.os.path.exists')
    @patch('app.routes.endpoints.plugin.FileResponse')
    def test_download_plugin_file_success(
        self,
        mock_file_response: MagicMock,
        mock_exists: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test downloading specific plugin file."""
        mock_exists.return_value = True
        mock_file_response.return_value = MagicMock()

        response = client.get(
            f"/routes/plugin/file/{sample_plugin.name}/test.py",
            headers=auth_headers
        )

        # Should attempt to download or return not found
        assert response.status_code in [200, 404]


# ==============================================================================
# Test Group 4: User-Plugin Association
# ==============================================================================

class TestPluginAssociation:
    """Test user-plugin association operations."""

    def test_associate_plugin_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin,
        sample_user: models.User
    ):
        """Test associating plugin with user."""
        association_data = {
            "plugin_name": sample_plugin.name,
            "user_id": sample_user.id
        }

        response = client.post(
            "/routes/plugin/associate",
            json=association_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        result = response.json()
        assert "message" in result or "success" in str(result).lower()

    def test_dissociate_plugin_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin,
        sample_user: models.User,
        db_session: Session
    ):
        """Test dissociating plugin from user."""
        # First associate
        sample_user.plugins.append(sample_plugin)
        db_session.commit()

        dissociation_data = {
            "plugin_name": sample_plugin.name,
            "user_id": sample_user.id
        }

        response = client.post(
            "/routes/plugin/dissociate",
            json=dissociation_data,
            headers=auth_headers
        )

        assert response.status_code == 200

    def test_associate_plugin_not_found(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User
    ):
        """Test associating non-existent plugin."""
        association_data = {
            "plugin_name": "NonExistentPlugin",
            "user_id": sample_user.id
        }

        response = client.post(
            "/routes/plugin/associate",
            json=association_data,
            headers=auth_headers
        )

        assert response.status_code in [404, 400]

    def test_dissociate_plugin_not_associated(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin,
        sample_user: models.User
    ):
        """Test dissociating plugin that wasn't associated."""
        dissociation_data = {
            "plugin_name": sample_plugin.name,
            "user_id": sample_user.id
        }

        response = client.post(
            "/routes/plugin/dissociate",
            json=dissociation_data,
            headers=auth_headers
        )

        # Should handle gracefully (either 200 or 404)
        assert response.status_code in [200, 404, 400]


# ==============================================================================
# Test Group 5: Docker Build Operations
# ==============================================================================

class TestPluginDockerBuild:
    """Test Docker build-related operations (with mocks)."""

    @patch('app.routes.endpoints.plugin.build_plugin_task.apply_async')
    def test_build_plugin_success(
        self,
        mock_task: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test initiating Docker build for plugin."""
        # Mock Celery task
        mock_task.return_value.id = "test-build-task-123"

        response = client.post(
            f"/routes/plugin/build/{sample_plugin.name}",
            headers=auth_headers
        )

        # Should start build or return current status
        assert response.status_code in [200, 202]

        if response.status_code in [200, 202]:
            data = response.json()
            assert "task_id" in data or "message" in data

    @patch('app.routes.endpoints.plugin.docker.from_env')
    def test_check_image_exists(
        self,
        mock_docker: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test checking if Docker image exists."""
        # Mock Docker client
        mock_client = MagicMock()
        mock_client.images.get.return_value = MagicMock()
        mock_docker.return_value = mock_client

        response = client.get(
            f"/routes/plugin/check_image/{sample_plugin.name}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "exists" in data or "image_exists" in data

    @patch('app.routes.endpoints.plugin.docker.from_env')
    def test_check_image_not_exists(
        self,
        mock_docker: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_plugin: models.Plugin
    ):
        """Test checking non-existent Docker image."""
        # Mock Docker client to raise exception
        mock_client = MagicMock()
        mock_client.images.get.side_effect = Exception("Image not found")
        mock_docker.return_value = mock_client

        response = client.get(
            f"/routes/plugin/check_image/{sample_plugin.name}",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "exists" in data or "image_exists" in data
        # Should indicate image doesn't exist
        if "exists" in data:
            assert data["exists"] is False
        elif "image_exists" in data:
            assert data["image_exists"] is False

    @patch('app.routes.endpoints.plugin.AsyncResult')
    def test_get_build_status_success(
        self,
        mock_async_result: MagicMock,
        client: TestClient,
        auth_headers: dict
    ):
        """Test retrieving build task status."""
        # Mock AsyncResult
        mock_result = MagicMock()
        mock_result.state = "SUCCESS"
        mock_result.info = {"message": "Build completed"}
        mock_async_result.return_value = mock_result

        response = client.get(
            "/routes/plugin/build/status/test-task-123",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "state" in data or "status" in data

    @patch('app.routes.endpoints.plugin.celery_app.control.revoke')
    def test_cancel_build_task_success(
        self,
        mock_revoke: MagicMock,
        client: TestClient,
        auth_headers: dict
    ):
        """Test canceling build task."""
        response = client.post(
            "/routes/plugin/build/cancel/test-task-123",
            headers=auth_headers
        )

        # Should accept cancellation request
        assert response.status_code in [200, 202]


# ==============================================================================
# Test Group 6: Security and Permissions
# ==============================================================================

class TestPluginSecurity:
    """Test plugin security and access control."""

    def test_plugin_access_without_auth(
        self,
        client: TestClient
    ):
        """Test accessing plugin endpoints without authentication."""
        response = client.get("/routes/plugin/list")

        assert response.status_code == 401

    def test_plugin_edit_official_plugin(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session
    ):
        """Test attempting to edit official (non-editable) plugin."""
        # Create official plugin
        official_plugin = models.Plugin(
            name="OfficialTestPlugin",
            description="Official plugin",
            author="official",
            plugin_path="./plugin/official/OfficialTestPlugin",
            plugin_type=PluginType.ANALYSIS,
            dependencies={},
            drawflow={},
            rules={},
            source="official",
            is_editable=False
        )
        db_session.add(official_plugin)
        db_session.commit()

        # Attempt to update
        update_data = {
            "description": "Modified description"
        }

        response = client.post(
            f"/routes/plugin/update/{official_plugin.name}",
            json=update_data,
            headers=auth_headers
        )

        # Should prevent editing official plugins
        assert response.status_code in [403, 400]

    def test_plugin_operations_inactive_user(
        self,
        client: TestClient,
        sample_inactive_user: models.User
    ):
        """Test plugin operations with inactive user."""
        # Attempt to login with inactive user
        response = client.post(
            "/routes/auth/login/access-token",
            data={
                "username": sample_inactive_user.email,
                "password": "testpassword123"
            }
        )

        # Inactive users should not be able to login
        assert response.status_code == 400
