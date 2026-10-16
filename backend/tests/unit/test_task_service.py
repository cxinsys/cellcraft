"""
Unit tests for the task service + SSE layer (app/task/service.py, app/task/sse.py),
extracted in PR-7.

Scope: pin the extracted business logic directly with the Celery / Docker /
filesystem boundaries mocked. The characterization/integration tests already
pin the HTTP contract; these tests give the new service functions and the SSE
streaming generator direct coverage.
"""
import pytest
from unittest.mock import patch, MagicMock

from fastapi import HTTPException

from app.core.exceptions import ValidationFailedError, ExternalServiceError
from app.task import service
from app.task import sse


def _user(username="tasksvc", user_id=11):
    u = MagicMock()
    u.id = user_id
    u.username = username
    return u


# ---------------------------------------------------------------------------
# SSE streaming: terminal states, timeout, cleanup
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
class TestStreamTaskStatus:
    """``sse.stream_task_status`` terminal/timeout/cleanup behavior."""

    async def _collect(self, task_id):
        return [frame async for frame in sse.stream_task_status(task_id)]

    @pytest.mark.parametrize("terminal", ["SUCCESS", "FAILURE", "REVOKED"])
    async def test_terminal_status_yields_once_and_stops(self, terminal):
        with patch("app.task.sse.get_task_info", return_value={"task_status": terminal}):
            frames = await self._collect("t1")
        assert frames == [terminal]

    async def test_missing_status_defaults_to_unknown_then_sleeps(self):
        # UNKNOWN is non-terminal -> yields once, then hits asyncio.sleep(5).
        # Patch sleep to raise so the loop terminates deterministically.
        with patch("app.task.sse.get_task_info", return_value={}), \
                patch("app.task.sse.asyncio.sleep", side_effect=StopAsyncIteration):
            frames = await self._collect("t2")
        # StopAsyncIteration is swallowed by the generator's except Exception.
        assert frames == ["UNKNOWN"]

    async def test_timeout_yields_timeout_frame(self):
        # Force the elapsed-time check to exceed the 3600s timeout immediately.
        times = iter([1000.0, 1000.0 + 4000])
        with patch("app.task.sse.get_task_info", return_value={"task_status": "RUNNING"}), \
                patch("app.task.sse.time.time", side_effect=lambda: next(times)):
            frames = await self._collect("t3")
        assert frames == ["TIMEOUT"]

    async def test_exception_is_swallowed(self):
        with patch("app.task.sse.get_task_info", side_effect=RuntimeError("boom")):
            frames = await self._collect("t4")
        assert frames == []


# ---------------------------------------------------------------------------
# get_task_status_simple
# ---------------------------------------------------------------------------

class TestGetTaskStatusSimple:
    def test_returns_status_shape(self):
        with patch("app.task.service.get_task_info", return_value={"task_status": "SUCCESS"}):
            result = service.get_task_status_simple(task_id="abc")
        assert result == {"task_id": "abc", "status": "SUCCESS"}

    def test_missing_status_defaults_unknown(self):
        with patch("app.task.service.get_task_info", return_value={}):
            result = service.get_task_status_simple(task_id="abc")
        assert result == {"task_id": "abc", "status": "UNKNOWN"}

    def test_blank_task_id_raises_400(self):
        with pytest.raises(ValidationFailedError) as exc:
            service.get_task_status_simple(task_id="   ")
        assert exc.value.status_code == 400
        assert exc.value.detail == "Invalid task_id"


# ---------------------------------------------------------------------------
# get_resources
# ---------------------------------------------------------------------------

class TestGetResources:
    def test_returns_status_dict(self):
        with patch("app.shared.resources.get_resource_status", return_value={"cpu": 4}):
            assert service.get_resources() == {"cpu": 4}

    def test_none_status_raises_503(self):
        with patch("app.shared.resources.get_resource_status", return_value=None):
            with pytest.raises(ExternalServiceError) as exc:
                service.get_resources()
        assert exc.value.status_code == 503
        assert exc.value.detail == "Resource monitoring unavailable"


# ---------------------------------------------------------------------------
# get_task_monitoring
# ---------------------------------------------------------------------------

class TestGetTaskMonitoring:
    def test_empty_returns_empty_list(self):
        db = MagicMock()
        with patch("app.task.service.crud_task.get_user_task_with_plugin", return_value=[]):
            assert service.get_task_monitoring(db=db, current_user=_user()) == []

    def test_filters_plugin_build_and_shapes_items(self):
        from datetime import datetime

        db = MagicMock()
        user = _user()

        def _task(task_type, tid):
            t = MagicMock()
            t.id = 1
            t.task_id = tid
            t.workflow_id = 5
            t.user_id = user.id
            t.status = "RUNNING"
            t.start_time = datetime(2026, 1, 1, 0, 0, 0)
            t.end_time = None
            t.task_type = task_type
            t.plugin_name = "P"
            t.plugin = None
            t.workflows = None
            return t

        tasks = [_task("compile", "c1"), _task("visualization", "v1"), _task("plugin_build", "b1")]
        with patch("app.task.service.crud_task.get_user_task_with_plugin", return_value=tasks):
            items = service.get_task_monitoring(db=db, current_user=user)

        assert len(items) == 2
        types = {i.task_type for i in items}
        assert "plugin_build" not in types
        assert types == {"compile", "visualization"}
        # compile task_title == plugin_name; visualization gets the -visualization suffix
        titles = {i.task_type: i.task_title for i in items}
        assert titles["compile"] == "P"
        assert titles["visualization"] == "P-visualization"


# ---------------------------------------------------------------------------
# delete_task
# ---------------------------------------------------------------------------

class TestDeleteTask:
    def test_returns_fixed_message(self):
        db = MagicMock()
        with patch("app.task.service.crud_task.delete_user_task", return_value=MagicMock()):
            result = service.delete_task(db=db, current_user=_user(), task_id="abc")
        assert result == {"message": "Task Deleted", "task_id": "abc"}

    def test_propagates_not_found(self):
        db = MagicMock()
        with patch(
            "app.task.service.crud_task.delete_user_task",
            side_effect=HTTPException(status_code=404, detail="Task abc not found"),
        ):
            with pytest.raises(HTTPException) as exc:
                service.delete_task(db=db, current_user=_user(), task_id="abc")
        assert exc.value.status_code == 404


# ---------------------------------------------------------------------------
# get_task_logs: auth + not-found
# ---------------------------------------------------------------------------

class TestGetTaskLogs:
    def test_task_not_found_raises_404(self):
        db = MagicMock()
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=None):
            with pytest.raises(HTTPException) as exc:
                service.get_task_logs(db=db, current_user=_user(), task_id="abc")
        assert exc.value.status_code == 404
        assert exc.value.detail == "Task not found"

    def test_other_users_task_raises_403(self):
        db = MagicMock()
        task = MagicMock()
        task.user_id = 999  # different from _user().id
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=task):
            with pytest.raises(HTTPException) as exc:
                service.get_task_logs(db=db, current_user=_user(), task_id="abc")
        assert exc.value.status_code == 403
        assert exc.value.detail == "Access denied"


# ---------------------------------------------------------------------------
# get_dag_structure: parser delegation + fallback + auth
# ---------------------------------------------------------------------------

class TestGetDagStructure:
    def _task(self, user):
        t = MagicMock()
        t.user_id = user.id
        t.task_id = "abc"
        t.workflow_id = 5
        t.algorithm_id = 1
        t.task_type = "compile"
        t.plugin_name = "P"
        t.status = "RUNNING"
        return t

    def test_native_parser_success(self):
        db = MagicMock()
        user = _user()
        dag = {"nodes": [{"id": "r1"}], "edges": [], "execution_sequence": ["r1"], "method": "native"}
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=self._task(user)), \
                patch("app.task.service.parse_snakefile_native", return_value=dag):
            result = service.get_dag_structure(db=db, current_user=user, task_id="abc")
        assert result["dag_structure"]["parsing_method"] == "native_native"
        assert result["task_info"]["task_id"] == "abc"

    def test_falls_back_to_legacy_parser(self):
        db = MagicMock()
        user = _user()
        legacy_parser = MagicMock()
        legacy_parser.parse_snakefile_with_logs.return_value = {
            "nodes": [{"id": "r1"}], "edges": [], "execution_sequence": ["r1"]
        }
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=self._task(user)), \
                patch("app.task.service.parse_snakefile_native", side_effect=Exception("native down")), \
                patch("app.task.service.SnakemakeDAGParser", return_value=legacy_parser):
            result = service.get_dag_structure(db=db, current_user=user, task_id="abc")
        assert result["dag_structure"]["parsing_method"] == "legacy"

    def test_task_not_found_raises_404(self):
        db = MagicMock()
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=None):
            with pytest.raises(HTTPException) as exc:
                service.get_dag_structure(db=db, current_user=_user(), task_id="abc")
        assert exc.value.status_code == 404

    def test_forbidden_raises_403(self):
        db = MagicMock()
        task = MagicMock()
        task.user_id = 999
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=task):
            with pytest.raises(HTTPException) as exc:
                service.get_dag_structure(db=db, current_user=_user(), task_id="abc")
        assert exc.value.status_code == 403


# ---------------------------------------------------------------------------
# get_rule_status: tracker success + fallback
# ---------------------------------------------------------------------------

class TestGetRuleStatus:
    def _task(self, user):
        t = MagicMock()
        t.user_id = user.id
        t.task_id = "abc"
        t.workflow_id = 5
        t.algorithm_id = 1
        t.task_type = "compile"
        t.plugin_name = "P"
        t.status = "RUNNING"
        return t

    def test_tracker_success(self):
        db = MagicMock()
        user = _user()
        dag = {"nodes": [{"id": "r1"}, {"id": "r2"}], "execution_sequence": ["r1", "r2"]}
        tracker = MagicMock()
        tracker.get_enhanced_progress_info.return_value = {
            "basic_progress": {
                "total_rules": 2, "completed_rules": 1, "failed_rules": 0,
                "running_rules": 1, "pending_rules": 0, "percentage": 50.0,
            },
            "rule_statuses": {"r1": "completed", "r2": "running"},
        }
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=self._task(user)), \
                patch("app.task.service.parse_snakefile_native", return_value=dag), \
                patch("app.task.service.SnakemakeRuleStatusTracker", return_value=tracker):
            result = service.get_rule_status(db=db, current_user=user, task_id="abc")
        assert result["rule_statuses"] == {"r1": "completed", "r2": "running"}
        assert result["progress"]["percentage"] == 50.0

    def test_actual_status_param_used(self):
        db = MagicMock()
        user = _user()
        dag = {"nodes": [{"id": "r1"}], "execution_sequence": ["r1"]}
        tracker = MagicMock()
        tracker.get_enhanced_progress_info.return_value = {
            "basic_progress": {"total_rules": 1, "completed_rules": 1, "percentage": 100.0},
            "rule_statuses": {"r1": "completed"},
        }
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=self._task(user)), \
                patch("app.task.service.parse_snakefile_native", return_value=dag), \
                patch("app.task.service.SnakemakeRuleStatusTracker", return_value=tracker) as mk:
            service.get_rule_status(db=db, current_user=user, task_id="abc", actual_status="COMPLETED")
        assert mk.call_args[1]["task_status"] == "COMPLETED"

    def test_tracker_failure_falls_back_to_pending(self):
        db = MagicMock()
        user = _user()
        dag = {"nodes": [{"id": "r1"}, {"id": "r2"}], "execution_sequence": ["r1", "r2"]}
        with patch("app.task.service.crud_task.get_task_by_task_id", return_value=self._task(user)), \
                patch("app.task.service.parse_snakefile_native", return_value=dag), \
                patch("app.task.service.SnakemakeRuleStatusTracker", side_effect=Exception("tracker down")):
            result = service.get_rule_status(db=db, current_user=user, task_id="abc")
        assert result["rule_statuses"] == {"r1": "pending", "r2": "pending"}
        assert result["progress"]["pending_rules"] == 2
        assert result["progress"]["completed_rules"] == 0
