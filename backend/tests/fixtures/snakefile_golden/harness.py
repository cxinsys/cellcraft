"""Run golden cases through the workflow service entry points (compile / visualize).

Only the database, Celery dispatch, cache and plugin lookup boundaries are mocked;
input extraction, validation and Snakefile generation run as in production.
Works against both the reference revision and the current code.
"""
import json
import os
from unittest.mock import MagicMock, patch

HERE = os.path.dirname(os.path.abspath(__file__))


def load_cases():
    with open(os.path.join(HERE, "cases.json"), encoding="utf-8") as f:
        return json.load(f)


def _user():
    user = MagicMock()
    user.id = 7
    user.username = "john"
    return user


def render(case, workdir):
    """Run ``case`` with ``workdir`` as the backend root and return the generated Snakefile."""
    from app.workflow import service
    from app.workflow.schemas import WorkflowCreate, WorkflowVisualizationRequest

    plugin_dir = os.path.join(HERE, "templates", case["plugin"])
    previous_cwd = os.getcwd()
    os.chdir(workdir)
    try:
        user_workflow = MagicMock()
        db = MagicMock()
        db.query.return_value.filter_by.return_value.first.return_value = None
        common = [
            patch.object(service.crud_workflow, "get_user_workflow", return_value=user_workflow),
            patch.object(service.crud_workflow, "update_workflow"),
            patch.object(service, "get_plugin_path", return_value=(plugin_dir, "official")),
            patch.object(service.process_data_task, "apply_async", return_value=MagicMock(id="tid")),
            patch.object(service, "get_task_info", return_value={"task_status": "PENDING"}),
        ]
        if case["kind"] == "analysis":
            user_workflow.workflow_info = case["workflow_info"]
            for p in common:
                p.start()
            try:
                service.compile_workflow(db=db, user=_user(), workflow=WorkflowCreate(
                    id=1, title="t", thumbnail=None, workflow_info=case["workflow_info"]))
            finally:
                patch.stopall()
            snakefile = "user/john/workflow_1/algorithm_2/Snakefile"
        else:
            results = os.path.join("user", "john", "workflow_1", "algorithm_2", "results")
            os.makedirs(results, exist_ok=True)
            for param in case["request"]["selectedVisualizationParams"]:
                if param.get("selectedFile"):
                    open(os.path.join(results, param["selectedFile"]), "w").close()
            user_workflow.workflow_info = {"drawflow": {"Home": {"data": {}}}}
            extra = [
                patch.object(service, "check_cache_with_expiry", return_value=(False, None)),
                patch.object(service, "maybe_cleanup_cache"),
            ]
            for p in common + extra:
                p.start()
            try:
                response = service.visualize_data(db=db, user=_user(),
                                                  workflow_request=WorkflowVisualizationRequest(**case["request"]))
            finally:
                patch.stopall()
            if not isinstance(response, dict):
                raise AssertionError(f"visualize_data failed: {getattr(response, 'body', response)}")
            snakefile = "user/john/workflow_1/visualization_7/Snakefile"
        with open(snakefile, encoding="utf-8") as f:
            return f.read()
    finally:
        os.chdir(previous_cwd)
