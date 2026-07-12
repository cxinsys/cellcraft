"""
Characterization tests for the workflow compile endpoint.

Purpose: pin the CURRENT behavior of ``POST /routes/workflow/compile``
(compileWorkflow, app/workflow/router.py lines 47-185) so the later
service-extraction refactor (PR-6) can prove behavior is unchanged.

IMPORTANT current-behavior note (preserved intentionally, do NOT "fix"):
The entire endpoint body is wrapped in a single ``try/except Exception`` that
re-raises as ``HTTPException(status_code=400, detail=str(e))`` (lines 181-185).
There is NO ``except HTTPException`` guard before it, so inner HTTPExceptions
(including the 404 "plugin not found" and the 400 "No algorithm nodes") are
caught by the outer handler and surfaced as HTTP 400. These tests pin that
observable 400 behavior rather than the status code named at the inner ``raise``.

Source of truth read while writing these tests:
- app/workflow/router.py :: compileWorkflow
- app/workflow/schemas.py :: WorkflowCreate

Pinned behavior:
- Success returns 200 with keys: message, task_ids, algorithm_ids,
  task_algorithm_mapping, results.
- A workflow with no algorithm nodes -> 400.
- A missing/unknown workflow id -> the ``if user_workflow:`` branch is skipped,
  so the function returns None -> FastAPI serializes 200 null (pinned below).
"""
import pytest
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app import models
from app.core.config import settings

COMPILE_URL = f"{settings.ROUTES_STR}/workflow/compile"


def _algorithm_workflow_info(plugin_name: str = "TestPlugin") -> dict:
    # extract_all_algorithms checks value.get("class") == "Algorithm" (capital A)
    # and reads selectedPlugin from value.get("data", {}) — not from value directly.
    return {
        "drawflow": {
            "Home": {
                "data": {
                    "1": {
                        "id": 1,
                        "name": "algorithm",
                        "class": "Algorithm",
                        "data": {
                            "selectedPlugin": {"name": plugin_name, "source": "local"},
                            "selectedPluginInputOutput": {},
                            "selectedPluginRules": {},
                            "files": [],
                        },
                    }
                }
            }
        }
    }


@pytest.mark.integration
@pytest.mark.characterization
class TestCharacterizationWorkflowCompile:
    """Freeze current compile behavior: success shape + error wrapping to 400."""

    @patch("app.workflow.router.process_data_task.apply_async")
    @patch("app.workflow.router.change_snakefile_parameter")
    @patch("app.workflow.router.get_plugin_path")
    @patch("app.workflow.router.get_task_info")
    def test_compile_success_response_shape(
        self,
        mock_get_task_info,
        mock_get_plugin_path,
        mock_change_snakefile,
        mock_apply_async,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session,
    ):
        """A single-algorithm workflow compiles and returns the current response keys."""
        workflow_info = _algorithm_workflow_info()
        workflow = models.Workflow(
            title="Char Compile Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id,
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        mock_get_plugin_path.return_value = ("./plugin/local/TestPlugin", "local")
        mock_change_snakefile.return_value = "/tmp/char_snakefile"
        mock_task = MagicMock()
        mock_task.id = "char-task-id-1"
        mock_apply_async.return_value = mock_task
        mock_get_task_info.return_value = {"task_status": "PENDING"}

        payload = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info,
        }
        response = client.post(COMPILE_URL, json=payload, headers=auth_headers)

        assert response.status_code == 200, response.text
        data = response.json()
        # Current top-level keys (pinned exactly).
        assert set(data.keys()) == {
            "message",
            "task_ids",
            "algorithm_ids",
            "task_algorithm_mapping",
            "results",
        }
        assert data["message"] == "Multiple tasks added to queue"
        assert data["task_ids"] == ["char-task-id-1"]
        assert data["algorithm_ids"] == [1]
        assert data["task_algorithm_mapping"] == {"char-task-id-1": 1}

    def test_compile_no_algorithm_nodes_returns_400(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session,
    ):
        """A workflow with no algorithm nodes surfaces as HTTP 400 (outer wrapper)."""
        workflow_info = {
            "drawflow": {"Home": {"data": {"1": {"id": 1, "name": "data", "class": "data"}}}}
        }
        workflow = models.Workflow(
            title="Char No-Algorithm Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id,
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        payload = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info,
        }
        response = client.post(COMPILE_URL, json=payload, headers=auth_headers)

        # The inner HTTPException(400, "No algorithm nodes found...") is caught by
        # the outer ``except Exception`` wrapper, which calls str(e) on the
        # HTTPException object. str(HTTPException) yields an empty string in
        # Starlette/FastAPI 0.78, so detail is "". Pin the actual observed shape.
        assert response.status_code == 400
        assert response.json()["detail"] == ""

    @patch("app.workflow.router.get_plugin_path")
    def test_compile_plugin_not_found_surfaces_as_400(
        self,
        mock_get_plugin_path,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User,
        db_session: Session,
    ):
        """Plugin resolution failure surfaces as HTTP 400 via the outer wrapper.

        The inner code raises HTTPException(404,...), but there is no
        ``except HTTPException`` guard, so the outer ``except Exception`` catches
        it and re-raises as 400. This is the CURRENT behavior being pinned.
        """
        workflow_info = _algorithm_workflow_info(plugin_name="NoSuchPlugin")
        workflow = models.Workflow(
            title="Char Bad-Plugin Workflow",
            workflow_info=workflow_info,
            user_id=sample_user.id,
        )
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        mock_get_plugin_path.side_effect = Exception("plugin missing")

        payload = {
            "id": workflow.id,
            "title": workflow.title,
            "thumbnail": None,
            "workflow_info": workflow_info,
        }
        response = client.post(COMPILE_URL, json=payload, headers=auth_headers)

        assert response.status_code == 400
        assert "detail" in response.json()

    def test_compile_unknown_workflow_returns_200_null(
        self,
        client: TestClient,
        auth_headers: dict,
    ):
        """When the workflow id is not owned/found, the ``if user_workflow`` branch
        is skipped and the function returns None -> FastAPI serializes 200 null.

        This pins a quirky-but-real current behavior (no explicit 404/400).
        """
        payload = {
            "id": 999999,
            "title": "ghost",
            "thumbnail": None,
            "workflow_info": {"drawflow": {"Home": {"data": {}}}},
        }
        response = client.post(COMPILE_URL, json=payload, headers=auth_headers)

        assert response.status_code == 200
        assert response.json() is None

    def test_compile_requires_authentication(
        self,
        client: TestClient,
    ):
        """Without auth headers the endpoint rejects with 401."""
        payload = {
            "id": 1,
            "title": "x",
            "thumbnail": None,
            "workflow_info": {"drawflow": {"Home": {"data": {}}}},
        }
        response = client.post(COMPILE_URL, json=payload)
        assert response.status_code == 401
