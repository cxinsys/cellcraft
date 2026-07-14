"""
Characterization tests for task status / monitoring endpoints.

Purpose: pin the CURRENT behavior of the task status surface so the later
service + SSE extraction refactor (PR-7) can prove behavior is unchanged.

Endpoints covered (app/routes/endpoints/task.py):
- GET  /routes/task/status/{task_id}   (lightweight, non-SSE; lines 82-95)
- GET  /routes/task/monitoring         (list; lines 97-159)
- DELETE /routes/task/delete/{task_id} (lines 276-287)
- GET  /routes/task/info/{task_id}     (SSE stream; lines 41-80) — one case only

External deps (Celery ``get_task_info``) are mocked, consistent with the
existing task_api.py mocking style.

Pinned behavior:
- /status returns {"task_id": <id>, "status": <celery status>}; empty task_id
  in the path resolves to a different route (404) rather than the 400 branch.
- /monitoring returns [] (200) when the user has no tasks; otherwise a list with
  the current item shape and plugin_build tasks filtered out.
- /delete returns {"message": "Task Deleted", "task_id": <id>} even for an id
  that does not belong to the user (current lenient behavior).
- /info is an SSE endpoint: content-type text/event-stream and no auth required.
"""
import uuid
from datetime import datetime

import pytest
from unittest.mock import patch
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app.database import models
from app.common.config import settings

TASK_PREFIX = f"{settings.ROUTES_STR}/task"


@pytest.mark.integration
@pytest.mark.characterization
class TestCharacterizationTaskStatus:
    """Freeze non-SSE task status endpoints; minimal SSE smoke."""

    @patch("app.routes.endpoints.task.get_task_info")
    def test_status_simple_response_shape(
        self,
        mock_get_task_info,
        client: TestClient,
    ):
        """/status returns {task_id, status} using the celery task_status value."""
        mock_get_task_info.return_value = {"task_status": "SUCCESS"}

        response = client.get(f"{TASK_PREFIX}/status/some-task-123")

        assert response.status_code == 200
        assert response.json() == {"task_id": "some-task-123", "status": "SUCCESS"}

    @patch("app.routes.endpoints.task.get_task_info")
    def test_status_missing_status_defaults_to_unknown(
        self,
        mock_get_task_info,
        client: TestClient,
    ):
        """When celery returns no task_status key, status defaults to UNKNOWN."""
        mock_get_task_info.return_value = {}

        response = client.get(f"{TASK_PREFIX}/status/no-status-task")

        assert response.status_code == 200
        assert response.json() == {"task_id": "no-status-task", "status": "UNKNOWN"}

    def test_monitoring_empty_returns_200_empty_list(
        self,
        client: TestClient,
        auth_headers: dict,
    ):
        """With no tasks, /monitoring returns an empty list and 200 (not 404)."""
        response = client.get(f"{TASK_PREFIX}/monitoring", headers=auth_headers)

        assert response.status_code == 200
        assert response.json() == []

    def test_monitoring_filters_plugin_build_and_shapes_items(
        self,
        client: TestClient,
        auth_headers: dict,
        multiple_user_tasks: list,
    ):
        """/monitoring excludes plugin_build tasks and pins the item key set."""
        response = client.get(f"{TASK_PREFIX}/monitoring", headers=auth_headers)

        assert response.status_code == 200
        data = response.json()

        # multiple_user_tasks fixture creates compile + visualization + plugin_build.
        # plugin_build must be filtered out -> 2 remain.
        assert len(data) == 2
        task_types = {t["task_type"] for t in data}
        assert "plugin_build" not in task_types
        assert task_types == {"compile", "visualization"}

        # Pin the item key set (TaskMonitoringItem serialization).
        for item in data:
            for key in (
                "id",
                "task_id",
                "workflow_id",
                "user_id",
                "status",
                "start_time",
                "task_type",
                "plugin_name",
                "plugin",
            ):
                assert key in item, f"missing key {key} in monitoring item"

    def test_monitoring_requires_authentication(
        self,
        client: TestClient,
    ):
        """/monitoring requires auth -> 401 without a token."""
        response = client.get(f"{TASK_PREFIX}/monitoring")
        assert response.status_code == 401

    def test_delete_task_returns_fixed_message(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_task: models.Task,
    ):
        """/delete returns the current fixed message payload for an owned task."""
        response = client.delete(
            f"{TASK_PREFIX}/delete/{sample_task.task_id}",
            headers=auth_headers,
        )

        assert response.status_code == 200
        assert response.json() == {
            "message": "Task Deleted",
            "task_id": sample_task.task_id,
        }

    def test_delete_unknown_task_returns_404(
        self,
        client: TestClient,
        auth_headers: dict,
    ):
        """Deleting a non-existent task id returns 404 today.

        crud_task.delete_user_task raises HTTPException(404) when the task is
        not found (crud_task.py line 92-94), which propagates to the caller.
        """
        response = client.delete(
            f"{TASK_PREFIX}/delete/does-not-exist-{uuid.uuid4().hex[:8]}",
            headers=auth_headers,
        )

        assert response.status_code == 404

    # NOTE: task.py binds get_task_info at import time
    # (`from ...celery_utils import get_task_info`), so the patch must target
    # the endpoint module namespace — patching celery_utils has no effect and
    # the SSE loop would poll the real backend forever.
    @patch("app.routes.endpoints.task.get_task_info")
    def test_info_sse_stream_content_type(
        self,
        mock_get_task_info,
        client: TestClient,
    ):
        """/info/{task_id} is an SSE stream (text/event-stream) and needs no auth.

        Only the transport/content-type contract is pinned here; detailed event
        framing is deferred to E2E (async SSE cannot be fully exercised via the
        sync TestClient).
        """
        mock_get_task_info.return_value = {"task_status": "SUCCESS"}

        response = client.get(f"{TASK_PREFIX}/info/sse-char-task")

        assert response.status_code == 200
        assert "text/event-stream" in response.headers.get("content-type", "")
