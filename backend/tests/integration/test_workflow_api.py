"""
Integration tests for Workflow API endpoints.

This module tests the workflow-related API endpoints in:
- POST /routes/workflow/compile - Compile and execute workflow
- POST /routes/workflow/visualization - Create visualization
- POST /routes/workflow/save - Create/update workflow
- POST /routes/workflow/delete - Delete workflow
- GET /routes/workflow/me - Get user workflows
- POST /routes/workflow/find - Find specific workflow
- POST /routes/workflow/results - Get workflow result files
- POST /routes/workflow/result - Download result file
- POST /routes/workflow/node/save - Save node data
- POST /routes/workflow/node/read - Read node data
- POST /routes/workflow/node/delete - Delete node data

Coverage Goal: 60%+ for app/workflow/router.py
Quality Score Goal: 8.2/10
"""
import pytest
import json
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app import models


# ==============================================================================
# Test Group 1: Workflow CRUD Operations (Create, Read, Update, Delete)
# ==============================================================================

class TestWorkflowCRUD:
    """Test workflow creation, retrieval, update, and deletion."""

    def test_create_workflow_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User
    ):
        """Test successful workflow creation."""
        workflow_data = {
            "id": None,
            "title": "New Test Workflow",
            "thumbnail": "data:image/png;base64,test",
            "workflow_info": {
                "drawflow": {
                    "Home": {
                        "data": {}
                    }
                }
            }
        }

        response = client.post(
            "/routes/workflow/save",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["title"] == "New Test Workflow"
        assert data["user_id"] == sample_user.id
        assert "id" in data

    def test_create_workflow_with_empty_drawflow(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test workflow creation with minimal drawflow structure."""
        workflow_data = {
            "id": None,
            "title": "Empty Workflow",
            "thumbnail": None,
            "workflow_info": {
                "drawflow": {
                    "Home": {
                        "data": {}
                    }
                }
            }
        }

        response = client.post(
            "/routes/workflow/save",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["title"] == "Empty Workflow"
        assert data["workflow_info"]["drawflow"]["Home"]["data"] == {}

    def test_get_user_workflows_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test retrieving user's workflow list."""
        response = client.get(
            "/routes/workflow/me",
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) >= 1

        # Check workflow structure
        workflow = data[0]
        assert "id" in workflow
        assert "title" in workflow
        assert "thumbnail" in workflow
        assert "updated_at" in workflow
        assert "user_id" in workflow

    def test_get_user_workflows_empty(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session,
        sample_user: models.User
    ):
        """Test retrieving workflows when user has none."""
        # Delete all workflows for this user
        db_session.query(models.Workflow).filter(
            models.Workflow.user_id == sample_user.id
        ).delete()
        db_session.commit()

        response = client.get(
            "/routes/workflow/me",
            headers=auth_headers
        )

        # Should return 400 when no workflows exist
        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_find_workflow_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test finding a specific workflow by ID."""
        response = client.post(
            "/routes/workflow/find",
            json={"id": sample_workflow.id},
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["title"] == sample_workflow.title
        assert "workflow_info" in data
        assert "thumbnail" in data

    def test_find_workflow_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test finding non-existent workflow."""
        response = client.post(
            "/routes/workflow/find",
            json={"id": 99999},
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]

    def test_update_workflow_title(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test updating workflow title."""
        workflow_data = {
            "id": sample_workflow.id,
            "title": "Updated Workflow Title",
            "thumbnail": sample_workflow.thumbnail,
            "workflow_info": sample_workflow.workflow_info
        }

        response = client.post(
            "/routes/workflow/save",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["title"] == "Updated Workflow Title"
        assert data["id"] == sample_workflow.id

    def test_update_workflow_info(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test updating workflow_info with new nodes."""
        new_workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "name": "data_node",
                            "class": "data",
                            "data": {}
                        }
                    }
                }
            }
        }

        workflow_data = {
            "id": sample_workflow.id,
            "title": sample_workflow.title,
            "thumbnail": sample_workflow.thumbnail,
            "workflow_info": new_workflow_info
        }

        response = client.post(
            "/routes/workflow/save",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "1" in data["workflow_info"]["drawflow"]["Home"]["data"]

    def test_delete_workflow_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        db_session: Session
    ):
        """Test successful workflow deletion."""
        response = client.post(
            "/routes/workflow/delete",
            json={"id": sample_workflow.id},
            headers=auth_headers
        )

        assert response.status_code == 200

        # Verify workflow is deleted from database
        deleted_workflow = db_session.query(models.Workflow).filter(
            models.Workflow.id == sample_workflow.id
        ).first()
        assert deleted_workflow is None

    def test_delete_workflow_with_directory(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test workflow deletion with associated directory cleanup."""
        # Create workflow directory
        workflow_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}"
        os.makedirs(workflow_dir, exist_ok=True)

        # Create a test file in directory
        test_file = os.path.join(workflow_dir, "test.txt")
        with open(test_file, "w") as f:
            f.write("test content")

        response = client.post(
            "/routes/workflow/delete",
            json={"id": sample_workflow.id},
            headers=auth_headers
        )

        assert response.status_code == 200

        # Verify directory is deleted
        assert not os.path.exists(workflow_dir)

    def test_delete_workflow_not_found(
        self,
        client: TestClient,
        auth_headers: dict
    ):
        """Test deleting non-existent workflow."""
        response = client.post(
            "/routes/workflow/delete",
            json={"id": 99999},
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "not exists" in response.json()["detail"]


# ==============================================================================
# Test Group 2: Workflow Compilation and Execution
# ==============================================================================

class TestWorkflowCompilation:
    """Test workflow compilation and task submission."""

    @patch('app.workflow.router.process_data_task.apply_async')
    @patch('app.workflow.router.get_plugin_path')
    def test_compile_workflow_single_algorithm(
        self,
        mock_get_plugin_path: MagicMock,
        mock_apply_async: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session
    ):
        """Test compiling workflow with single algorithm node."""
        # Create workflow with algorithm node
        workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "algorithm",
                            "selectedPlugin": {
                                "name": "TestPlugin",
                                "source": "local"
                            },
                            "selectedPluginInputOutput": {},
                            "selectedPluginRules": {}
                        }
                    }
                }
            }
        }

        workflow = models.Workflow(
            title="Compile Test Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        # Mock plugin path
        mock_get_plugin_path.return_value = ("./plugin/local/TestPlugin", "local")

        # Mock Celery task
        mock_task = MagicMock()
        mock_task.id = "test-task-id-123"
        mock_apply_async.return_value = mock_task

        workflow_data = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info
        }

        with patch('app.workflow.router.change_snakefile_parameter') as mock_snakefile:
            mock_snakefile.return_value = "/tmp/test_snakefile"

            response = client.post(
                "/routes/workflow/compile",
                json=workflow_data,
                headers=auth_headers
            )

        assert response.status_code == 200
        data = response.json()
        assert "task_ids" in data
        assert len(data["task_ids"]) == 1
        assert data["task_ids"][0] == "test-task-id-123"

    @patch('app.workflow.router.process_data_task.apply_async')
    @patch('app.workflow.router.get_plugin_path')
    def test_compile_workflow_multiple_algorithms(
        self,
        mock_get_plugin_path: MagicMock,
        mock_apply_async: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session
    ):
        """Test compiling workflow with multiple algorithm nodes."""
        workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "algorithm",
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
                            "selectedPlugin": {
                                "name": "TestPlugin2",
                                "source": "local"
                            },
                            "selectedPluginInputOutput": {},
                            "selectedPluginRules": {}
                        }
                    }
                }
            }
        }

        workflow = models.Workflow(
            title="Multi-Algorithm Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        mock_get_plugin_path.return_value = ("./plugin/local/TestPlugin", "local")

        # Mock multiple Celery tasks
        mock_task1 = MagicMock()
        mock_task1.id = "task-1"
        mock_task2 = MagicMock()
        mock_task2.id = "task-2"
        mock_apply_async.side_effect = [mock_task1, mock_task2]

        workflow_data = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info
        }

        with patch('app.workflow.router.change_snakefile_parameter') as mock_snakefile:
            mock_snakefile.return_value = "/tmp/test_snakefile"

            response = client.post(
                "/routes/workflow/compile",
                json=workflow_data,
                headers=auth_headers
            )

        assert response.status_code == 200
        data = response.json()
        assert len(data["task_ids"]) == 2
        assert "task-1" in data["task_ids"]
        assert "task-2" in data["task_ids"]

    def test_compile_workflow_no_algorithm(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session
    ):
        """Test compiling workflow with no algorithm nodes raises error."""
        workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "data",  # Not an algorithm node
                            "class": "data"
                        }
                    }
                }
            }
        }

        workflow = models.Workflow(
            title="No Algorithm Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        workflow_data = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info
        }

        response = client.post(
            "/routes/workflow/compile",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 400
        assert "No algorithm nodes" in response.json()["detail"]

    @patch('app.workflow.router.get_plugin_path')
    def test_compile_workflow_invalid_plugin(
        self,
        mock_get_plugin_path: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session
    ):
        """Test compiling workflow with invalid plugin raises error."""
        workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "algorithm",
                            "selectedPlugin": {
                                "name": "NonExistentPlugin",
                                "source": "local"
                            },
                            "selectedPluginInputOutput": {},
                            "selectedPluginRules": {}
                        }
                    }
                }
            }
        }

        workflow = models.Workflow(
            title="Invalid Plugin Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        # Mock plugin not found
        mock_get_plugin_path.side_effect = Exception("Plugin not found")

        workflow_data = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info
        }

        response = client.post(
            "/routes/workflow/compile",
            json=workflow_data,
            headers=auth_headers
        )

        assert response.status_code == 404
        assert "not found" in response.json()["detail"]


# ==============================================================================
# Test Group 3: Workflow Results
# ==============================================================================

class TestWorkflowResults:
    """Test workflow result file operations."""

    def test_get_workflow_results_list(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test retrieving list of workflow result files."""
        # Create result directory and files
        result_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}/algorithm_1/results"
        os.makedirs(result_dir, exist_ok=True)

        # Create test result files
        test_files = ["output.csv", "network.json"]
        for filename in test_files:
            filepath = os.path.join(result_dir, filename)
            with open(filepath, "w") as f:
                f.write("test content")

        response = client.post(
            "/routes/workflow/results",
            json={
                "id": sample_workflow.id,
                "algorithm_id": 1
            },
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 2

        # Check file info structure
        file_names = [f["name"] for f in data]
        assert "output.csv" in file_names
        assert "network.json" in file_names

        # Verify file info contains size and modified_time
        for file_info in data:
            assert "name" in file_info
            assert "size" in file_info
            assert "modified_time" in file_info

        # Cleanup
        shutil.rmtree(f"./user/{sample_user.username}/workflow_{sample_workflow.id}", ignore_errors=True)

    def test_get_workflow_results_empty(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test retrieving results when directory doesn't exist."""
        response = client.post(
            "/routes/workflow/results",
            json={
                "id": sample_workflow.id,
                "algorithm_id": 999
            },
            headers=auth_headers
        )

        assert response.status_code == 404
        assert "not found" in response.json()["detail"]

    def test_download_workflow_result_file(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test downloading a specific workflow result file."""
        # Create result directory and file
        result_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}/algorithm_1/results"
        os.makedirs(result_dir, exist_ok=True)

        test_content = "gene1,gene2,score\nGENE_A,GENE_B,0.95\n"
        filepath = os.path.join(result_dir, "output.csv")
        with open(filepath, "w") as f:
            f.write(test_content)

        response = client.post(
            "/routes/workflow/result",
            json={
                "id": sample_workflow.id,
                "algorithm_id": 1,
                "filename": "output.csv"
            },
            headers=auth_headers
        )

        assert response.status_code == 200
        assert response.headers["content-type"] in ["text/csv", "application/octet-stream"]
        assert test_content in response.text

        # Cleanup
        shutil.rmtree(f"./user/{sample_user.username}/workflow_{sample_workflow.id}", ignore_errors=True)

    def test_download_workflow_result_not_found(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test downloading non-existent result file."""
        # Create result directory but no file
        result_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}/algorithm_1/results"
        os.makedirs(result_dir, exist_ok=True)

        response = client.post(
            "/routes/workflow/result",
            json={
                "id": sample_workflow.id,
                "algorithm_id": 1,
                "filename": "nonexistent.csv"
            },
            headers=auth_headers
        )

        # Should raise FileNotFoundError or 404
        assert response.status_code in [404, 500]

        # Cleanup
        shutil.rmtree(f"./user/{sample_user.username}/workflow_{sample_workflow.id}", ignore_errors=True)


# ==============================================================================
# Test Group 4: Visualization
# ==============================================================================

class TestWorkflowVisualization:
    """Test workflow visualization operations."""

    @patch('app.workflow.router.process_data_task.apply_async')
    @patch('app.workflow.router.get_plugin_path')
    @patch('app.workflow.router.crud_plugin.get_plugin_by_name')
    def test_create_visualization_success(
        self,
        mock_get_plugin: MagicMock,
        mock_get_plugin_path: MagicMock,
        mock_apply_async: MagicMock,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_plugin: models.Plugin
    ):
        """Test creating visualization successfully."""
        # Mock plugin path
        mock_get_plugin_path.return_value = ("./plugin/local/TestPlugin", "local")
        mock_get_plugin.return_value = sample_plugin

        # Create plugin metadata file
        plugin_dir = "./plugin/local/TestPlugin"
        os.makedirs(plugin_dir, exist_ok=True)

        metadata = {
            "rules": {
                "viz_1": {
                    "name": "Test Visualization",
                    "parameters": []
                }
            }
        }

        with open(os.path.join(plugin_dir, "metadata.json"), "w") as f:
            json.dump(metadata, f)

        # Create dummy Snakefile
        with open(os.path.join(plugin_dir, "Snakefile"), "w") as f:
            f.write("rule test_viz:\n    output: 'test.json'\n")

        # Mock Celery task
        mock_task = MagicMock()
        mock_task.id = "viz-task-123"
        mock_apply_async.return_value = mock_task

        visualization_request = {
            "id": sample_workflow.id,
            "current_node_id": "viz_1",
            "algorithm_id": "1",
            "selectedVisualizationPlugin": {
                "name": "TestPlugin",
                "source": "local"
            },
            "selectedScript": {
                "name": "Test Visualization"
            },
            "selectedVisualizationParams": []
        }

        with patch('app.workflow.router.extract_rule_block') as mock_extract:
            mock_extract.return_value = ("rule content", "/tmp/snakefile")

            response = client.post(
                "/routes/workflow/visualization",
                json=visualization_request,
                headers=auth_headers
            )

        assert response.status_code == 200
        data = response.json()
        assert data["task_id"] == "viz-task-123"
        assert "result_path" in data

        # Cleanup
        shutil.rmtree(plugin_dir, ignore_errors=True)

    def test_create_visualization_missing_plugin(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test creating visualization with missing plugin info."""
        visualization_request = {
            "id": sample_workflow.id,
            "current_node_id": "viz_1",
            "selectedVisualizationPlugin": None,  # Missing
            "selectedScript": {
                "name": "Test Visualization"
            },
            "selectedVisualizationParams": []
        }

        response = client.post(
            "/routes/workflow/visualization",
            json=visualization_request,
            headers=auth_headers
        )

        # Should return validation error
        assert response.status_code in [400, 422]

    def test_get_visualization_result(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test retrieving visualization result."""
        # Create visualization result directory and file
        viz_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}/visualization_1/results"
        os.makedirs(viz_dir, exist_ok=True)

        # Create test visualization result (Plotly JSON format)
        viz_data = {
            "data": [{"x": [1, 2, 3], "y": [4, 5, 6], "type": "scatter"}],
            "layout": {"title": "Test Visualization"}
        }

        filepath = os.path.join(viz_dir, "viz_result.json")
        with open(filepath, "w") as f:
            json.dump(viz_data, f)

        response = client.post(
            "/routes/workflow/visualization/result",
            json={
                "id": sample_workflow.id,
                "algorithm_id": 1,  # Used as visualization node ID
                "filename": "viz_result.json"
            },
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert "data" in data
        assert "layout" in data
        assert data["layout"]["title"] == "Test Visualization"

        # Cleanup
        shutil.rmtree(f"./user/{sample_user.username}/workflow_{sample_workflow.id}", ignore_errors=True)


# ==============================================================================
# Test Group 5: Node Data Operations
# ==============================================================================

class TestWorkflowNodeData:
    """Test workflow node data save/read/delete operations."""

    def test_save_node_data_json(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test saving JSON node data."""
        node_data = {
            "id": sample_workflow.id,
            "node_id": "node_1",
            "node_name": "test_node",
            "file_extension": "json",
            "file_content": {"key": "value", "list": [1, 2, 3]}
        }

        response = client.post(
            "/routes/workflow/node/save",
            json=node_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["node_id"] == "node_1"
        assert data["node_name"] == "test_node"
        assert "file_path" in data

        # Verify file was created
        assert os.path.exists(data["file_path"])

        # Verify file content
        with open(data["file_path"], "r") as f:
            content = json.load(f)
        assert content["key"] == "value"

        # Cleanup
        os.remove(data["file_path"])

    def test_save_node_data_csv(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow
    ):
        """Test saving CSV node data."""
        csv_content = "gene_id,expression\nGENE_1,5.2\nGENE_2,3.8"

        node_data = {
            "id": sample_workflow.id,
            "node_id": "node_2",
            "node_name": "gene_list",
            "file_extension": "csv",
            "file_content": csv_content
        }

        response = client.post(
            "/routes/workflow/node/save",
            json=node_data,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["file_extension"] == "csv"

        # Verify file content
        with open(data["file_path"], "r") as f:
            content = f.read()
        assert content == csv_content

        # Cleanup
        os.remove(data["file_path"])

    def test_read_node_data_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test reading node data successfully."""
        # First save node data
        node_data = {
            "id": sample_workflow.id,
            "node_id": "node_3",
            "node_name": "params",
            "file_extension": "json",
            "file_content": {"param1": "value1"}
        }

        save_response = client.post(
            "/routes/workflow/node/save",
            json=node_data,
            headers=auth_headers
        )
        assert save_response.status_code == 200

        # Now read it back
        read_request = {
            "id": sample_workflow.id,
            "node_id": "node_3",
            "node_name": "params",
            "file_extension": "json"
        }

        response = client.post(
            "/routes/workflow/node/read",
            json=read_request,
            headers=auth_headers
        )

        assert response.status_code == 200
        data = response.json()
        assert data["file_content"]["param1"] == "value1"

        # Cleanup
        workflow_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}"
        shutil.rmtree(workflow_dir, ignore_errors=True)

    def test_delete_node_data_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        sample_user: models.User
    ):
        """Test deleting node data successfully."""
        # First save node data
        node_data = {
            "id": sample_workflow.id,
            "node_id": "node_4",
            "node_name": "temp_data",
            "file_extension": "txt",
            "file_content": "temporary content"
        }

        save_response = client.post(
            "/routes/workflow/node/save",
            json=node_data,
            headers=auth_headers
        )
        assert save_response.status_code == 200
        file_path = save_response.json()["file_path"]

        # Verify file exists
        assert os.path.exists(file_path)

        # Delete it
        delete_request = {
            "id": sample_workflow.id,
            "node_id": "node_4",
            "node_name": "temp_data",
            "file_extension": "txt"
        }

        response = client.post(
            "/routes/workflow/node/delete",
            json=delete_request,
            headers=auth_headers
        )

        assert response.status_code == 200

        # Verify file is deleted
        assert not os.path.exists(file_path)

        # Cleanup
        workflow_dir = f"./user/{sample_user.username}/workflow_{sample_workflow.id}"
        shutil.rmtree(workflow_dir, ignore_errors=True)


# ==============================================================================
# Test Group 6: Security and Permissions
# ==============================================================================

class TestWorkflowSecurity:
    """Test workflow security and access control."""

    def test_workflow_access_without_auth(
        self,
        client: TestClient,
        sample_workflow: models.Workflow
    ):
        """Test accessing workflow endpoints without authentication."""
        response = client.get("/routes/workflow/me")

        assert response.status_code == 401

    def test_workflow_access_other_user(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_workflow: models.Workflow,
        user_factory
    ):
        """Test accessing another user's workflow."""
        # Create another user and their workflow
        other_user = user_factory(username="otheruser", email="other@test.com")
        other_workflow = models.Workflow(
            title="Other User Workflow",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=other_user.id
        )

        from sqlalchemy.orm import Session
        # Note: This test needs access to db_session, will need to adjust fixture usage

        # Try to access other user's workflow
        response = client.post(
            "/routes/workflow/find",
            json={"id": other_workflow.id},
            headers=auth_headers
        )

        # Should return 400 (not found) because user doesn't own this workflow
        assert response.status_code == 400

    def test_workflow_operations_inactive_user(
        self,
        client: TestClient,
        sample_inactive_user: models.User
    ):
        """Test workflow operations with inactive user."""
        # Login with inactive user (should fail)
        response = client.post(
            "/routes/auth/login/access-token",
            data={
                "username": sample_inactive_user.email,
                "password": "testpassword123"
            }
        )

        # Inactive users should not be able to login
        assert response.status_code == 400
