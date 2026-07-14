"""
Unit tests for the workflow service layer (app/workflow/service.py) extracted in
PR-6.

Scope: pin the extracted business logic directly, with the Celery dispatch,
filesystem, and cache dependencies mocked (using ``tmp_path`` for the real
directory operations the compile pipeline performs). The real workflow parsers
(``extract_all_algorithms``, ``extract_file_sources`` etc.) run unmocked; only
the compiler Snakefile generator and the cache/Celery boundaries are stubbed.

The characterization/integration tests already pin the HTTP contract; these
tests give the new service functions and the compile-pipeline step functions
(``_resolve_plugin_paths``, ``_prepare_task_workspace``,
``_generate_compile_snakefile``, ``_resolve_compile_resources``,
``_dispatch_compile_task``) direct coverage — including the visualization cache
hit and cache miss paths.
"""
import os
import json
import pytest
from unittest.mock import patch, MagicMock

from fastapi import HTTPException
from fastapi.responses import JSONResponse

from app.workflow import service
from app.workflow.schemas import (
    WorkflowCreate, WorkflowResult, WorkflowVisualizationRequest,
    WorkflowNodeFileCreate, WorkflowNodeFileRead, WorkflowNodeFileDelete,
)


def _user(username="svcuser", user_id=7):
    u = MagicMock()
    u.id = user_id
    u.username = username
    return u


def _algorithm_workflow_info(plugin_name="TestPlugin", source="local"):
    return {
        "drawflow": {
            "Home": {
                "data": {
                    "1": {
                        "id": 1,
                        "name": "algorithm",
                        "class": "Algorithm",
                        "data": {
                            "selectedPlugin": {"name": plugin_name, "source": source},
                            "selectedPluginInputOutput": {},
                            "selectedPluginRules": {},
                            "files": [],
                        },
                    }
                }
            }
        }
    }


# ---------------------------------------------------------------------------
# Compile pipeline step functions (direct coverage, tmp_path where real IO)
# ---------------------------------------------------------------------------

class TestCompileStepFunctions:
    """The private pipeline steps that compile_workflow orchestrates."""

    def test_resolve_plugin_paths_uses_frontend_source(self):
        algorithm = {"selectedPlugin": {"name": "P", "source": "official"}}
        with patch("app.workflow.service.get_plugin_path") as mock_gpp:
            mock_gpp.return_value = ("/plug/official/P", "official")
            plugin_path, dep_path = service._resolve_plugin_paths(algorithm=algorithm)

        mock_gpp.assert_called_once_with("P", "official")
        assert plugin_path == "/plug/official/P"
        assert dep_path == os.path.join("/plug/official/P", "dependency")

    def test_resolve_plugin_paths_autodetects_without_source(self):
        algorithm = {"selectedPlugin": {"name": "P"}}
        with patch("app.workflow.service.get_plugin_path") as mock_gpp:
            mock_gpp.return_value = ("/plug/local/P", "local")
            service._resolve_plugin_paths(algorithm=algorithm)
        # No source -> single-arg auto-detect call
        mock_gpp.assert_called_once_with("P")

    def test_resolve_plugin_paths_missing_plugin_raises_404(self):
        algorithm = {"selectedPlugin": {"name": "Ghost", "source": "local"}}
        with patch("app.workflow.service.get_plugin_path", side_effect=Exception("nope")):
            with pytest.raises(HTTPException) as exc:
                service._resolve_plugin_paths(algorithm=algorithm)
        assert exc.value.status_code == 404
        assert "Plugin 'Ghost' not found" in exc.value.detail

    def test_prepare_task_workspace_creates_missing_dir(self, tmp_path):
        task_path = str(tmp_path / "workflow_1" / "algorithm_1")
        user_workflow = MagicMock()
        service._prepare_task_workspace(
            user_workflow=user_workflow,
            user_path=str(tmp_path) + "/",
            workflow_id=1,
            algorithm={"id": 1},
            user_workflow_task_path=task_path,
        )
        assert os.path.isdir(task_path)

    def test_prepare_task_workspace_cleans_existing_results(self, tmp_path):
        task_path = str(tmp_path / "workflow_1" / "algorithm_1")
        os.makedirs(task_path)
        user_workflow = MagicMock()
        user_workflow.workflow_info = _algorithm_workflow_info()
        with patch("app.workflow.service.cleanup_task_results") as mock_cleanup, \
             patch("app.workflow.service.find_connected_visualization_nodes", return_value=[]):
            mock_cleanup.return_value = {"success": True, "files_removed": ["a"], "symlinks_removed": []}
            service._prepare_task_workspace(
                user_workflow=user_workflow,
                user_path=str(tmp_path) + "/",
                workflow_id=1,
                algorithm={"id": 1},
                user_workflow_task_path=task_path,
            )
        # Existing dir -> cleanup path taken (not makedirs)
        mock_cleanup.assert_called_once()

    def test_generate_compile_snakefile_delegates_to_compiler(self):
        with patch("app.workflow.service.change_snakefile_parameter") as mock_change:
            mock_change.return_value = "/tmp/task/Snakefile"
            out = service._generate_compile_snakefile(
                plugin_path="/plug/P",
                user_workflow_task_path="/tmp/task",
                user_input={"user_name": "x"},
                plugin_params={},
            )
        assert out == "/tmp/task/Snakefile"
        mock_change.assert_called_once_with(
            os.path.join("/plug/P", "Snakefile"),
            "/tmp/task/Snakefile",
            {"user_name": "x"},
            {},
        )

    def test_resolve_compile_resources_gpu_from_device_count(self):
        db = MagicMock()
        plugin = MagicMock()
        plugin.use_gpu = True
        db.query.return_value.filter_by.return_value.first.return_value = plugin
        rtype, slots = service._resolve_compile_resources(
            db=db, selected_plugin_name="P", plugin_params={"number of devices": "3"}
        )
        assert rtype == "gpu"
        assert slots == 3

    def test_resolve_compile_resources_defaults_cpu_when_no_plugin(self):
        db = MagicMock()
        db.query.return_value.filter_by.return_value.first.return_value = None
        with patch("app.core.config.settings") as mock_settings:
            mock_settings.RESOURCE_DEFAULT_CPU_SLOTS = 2
            mock_settings.RESOURCE_DEFAULT_GPU_SLOTS = 8
            rtype, slots = service._resolve_compile_resources(
                db=db, selected_plugin_name="P", plugin_params={}
            )
        assert rtype == "cpu"
        assert slots == 2

    def test_dispatch_compile_task_preserves_wire_kwargs(self):
        user = _user()
        with patch("app.workflow.service.process_data_task.apply_async") as mock_async:
            mock_async.return_value = MagicMock(id="task-1")
            task_id = service._dispatch_compile_task(
                user=user,
                workflow_id=5,
                algorithm={"id": 9},
                user_snakefile_path="/tmp/Snakefile",
                selected_plugin_name="P",
                target_list=["out.csv"],
                resource_type="cpu",
                resource_slots=1,
            )
        assert task_id == "task-1"
        _, kwargs = mock_async.call_args
        assert kwargs["args"] == [user.username, "/tmp/Snakefile", "P", ["out.csv"]]
        assert kwargs["kwargs"] == {
            "user_id": user.id,
            "workflow_id": 5,
            "algorithm_id": 9,
            "plugin_name": "P",
            "task_type": "compile",
            "resource_type": "cpu",
            "resource_slots": 1,
        }
        assert kwargs["ignore_result"] is False


# ---------------------------------------------------------------------------
# compile_workflow: success + error wrapping
# ---------------------------------------------------------------------------

class TestCompileWorkflow:
    """Full compile orchestration with real parsers, mocked IO/Celery boundaries."""

    def test_compile_success_shape(self, tmp_path):
        user = _user()
        db = MagicMock()
        wf_info = _algorithm_workflow_info()
        user_workflow = MagicMock()
        user_workflow.workflow_info = wf_info
        db.query.return_value.filter_by.return_value.first.return_value = None  # no plugin -> cpu default

        workflow = WorkflowCreate(id=1, title="t", thumbnail=None, workflow_info=wf_info)

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.crud_workflow.update_workflow"), \
             patch("app.workflow.service.get_plugin_path", return_value=("/plug/local/TestPlugin", "local")), \
             patch("app.workflow.service.generate_user_input", return_value={}), \
             patch("app.workflow.service.generate_plugin_params", return_value={}), \
             patch("app.workflow.service.extract_target_data", return_value=["out.csv"]), \
             patch("app.workflow.service.extract_file_sources", return_value={}), \
             patch("app.workflow.service.os.path.exists", return_value=False), \
             patch("app.workflow.service.os.makedirs"), \
             patch("app.workflow.service.change_snakefile_parameter", return_value="/tmp/Snakefile"), \
             patch("app.workflow.service.process_data_task.apply_async", return_value=MagicMock(id="tid-1")), \
             patch("app.workflow.service.get_task_info", return_value={"task_status": "PENDING"}), \
             patch("app.core.config.settings") as mock_settings:
            mock_settings.RESOURCE_DEFAULT_CPU_SLOTS = 1
            mock_settings.RESOURCE_DEFAULT_GPU_SLOTS = 8
            result = service.compile_workflow(db=db, user=user, workflow=workflow)

        assert set(result.keys()) == {
            "message", "task_ids", "algorithm_ids",
            "task_algorithm_mapping", "results",
        }
        assert result["message"] == "Multiple tasks added to queue"
        assert result["task_ids"] == ["tid-1"]
        assert result["algorithm_ids"] == [1]
        assert result["task_algorithm_mapping"] == {"tid-1": 1}

    def test_compile_no_algorithm_wraps_to_400(self):
        user = _user()
        db = MagicMock()
        wf_info = {"drawflow": {"Home": {"data": {"1": {"id": 1, "class": "data"}}}}}
        user_workflow = MagicMock()
        user_workflow.workflow_info = wf_info
        workflow = WorkflowCreate(id=1, title="t", thumbnail=None, workflow_info=wf_info)

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.crud_workflow.update_workflow"):
            with pytest.raises(HTTPException) as exc:
                service.compile_workflow(db=db, user=user, workflow=workflow)
        # Inner 400 is caught by the outer wrapper; str(HTTPException) is "".
        assert exc.value.status_code == 400

    def test_compile_plugin_not_found_wraps_to_400(self):
        user = _user()
        db = MagicMock()
        wf_info = _algorithm_workflow_info(plugin_name="Ghost")
        user_workflow = MagicMock()
        user_workflow.workflow_info = wf_info
        workflow = WorkflowCreate(id=1, title="t", thumbnail=None, workflow_info=wf_info)

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.crud_workflow.update_workflow"), \
             patch("app.workflow.service.get_plugin_path", side_effect=Exception("missing")):
            with pytest.raises(HTTPException) as exc:
                service.compile_workflow(db=db, user=user, workflow=workflow)
        # Inner 404 is re-wrapped as 400 by the outer except (pinned quirk).
        assert exc.value.status_code == 400

    def test_compile_unknown_workflow_returns_none(self):
        user = _user()
        db = MagicMock()
        workflow = WorkflowCreate(id=999, title="ghost", thumbnail=None,
                                  workflow_info={"drawflow": {"Home": {"data": {}}}})
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=None):
            result = service.compile_workflow(db=db, user=user, workflow=workflow)
        assert result is None


# ---------------------------------------------------------------------------
# visualize_data: cache hit + cache miss
# ---------------------------------------------------------------------------

def _viz_request(workflow_id=1):
    return WorkflowVisualizationRequest(
        id=workflow_id,
        current_node_id="viz_1",
        algorithm_id="1",
        selectedVisualizationPlugin={"name": "TestPlugin", "source": "local"},
        selectedScript={"name": "Test Visualization"},
        selectedVisualizationParams=[],
    )


class TestVisualizeData:
    """Visualization pipeline: cache-hit short circuit vs full dispatch."""

    def test_visualization_cache_hit_returns_cached(self, tmp_path):
        user = _user()
        db = MagicMock()
        user_workflow = MagicMock()
        user_workflow.workflow_info = _algorithm_workflow_info()

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.get_plugin_path", return_value=("/plug/local/TestPlugin", "local")), \
             patch("app.workflow.service.generate_cache_key", return_value="ck-1"), \
             patch("app.workflow.service.maybe_cleanup_cache"), \
             patch("app.workflow.service.check_cache_with_expiry", return_value=(True, "/cache/ck-1.json")), \
             patch("app.workflow.service.os.makedirs"), \
             patch("app.workflow.service.create_symbolic_link", return_value=True), \
             patch("app.workflow.service.update_cache_link_location"), \
             patch("app.workflow.service.process_data_task.apply_async") as mock_async:
            result = service.visualize_data(db=db, user=user, workflow_request=_viz_request())

        assert result == {
            "message": "Visualization result from cache",
            "result_path": "visualization_result.json",
            "cached": True,
        }
        # Cache hit must NOT dispatch a Celery task.
        mock_async.assert_not_called()

    def test_visualization_cache_miss_dispatches_task(self, tmp_path):
        user = _user()
        db = MagicMock()
        user_workflow = MagicMock()
        user_workflow.workflow_info = _algorithm_workflow_info()

        plugin_info = MagicMock()
        metadata = {"rules": {"r1": {"name": "Test Visualization"}}}

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.get_plugin_path", return_value=("/plug/local/TestPlugin", "local")), \
             patch("app.workflow.service.generate_cache_key", return_value="ck-2"), \
             patch("app.workflow.service.maybe_cleanup_cache"), \
             patch("app.workflow.service.check_cache_with_expiry", return_value=(False, None)), \
             patch("app.workflow.service.os.makedirs"), \
             patch("app.workflow.service.crud_plugin.get_plugin_by_name", return_value=plugin_info), \
             patch("builtins.open", new_callable=MagicMock), \
             patch("app.workflow.service.json.load", return_value=metadata), \
             patch("app.workflow.service.extract_rule_block", return_value=("content", "/tmp/viz/Snakefile")), \
             patch("app.workflow.service.process_data_task.apply_async", return_value=MagicMock(id="viz-tid")), \
             patch("app.workflow.service.get_task_info", return_value={"task_status": "PENDING"}):
            result = service.visualize_data(db=db, user=user, workflow_request=_viz_request())

        assert result["message"] == "Visualization task added to queue"
        assert result["task_id"] == "viz-tid"
        assert result["result_path"] == "visualization_result.json"
        assert result["cached"] is False

    def test_visualization_missing_plugin_raises_validation(self):
        user = _user()
        db = MagicMock()
        user_workflow = MagicMock()
        req = _viz_request()
        req.selectedVisualizationPlugin = None  # trigger ValidationError branch

        from app.core.exceptions import ValidationError
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow):
            with pytest.raises(ValidationError):
                service.visualize_data(db=db, user=user, workflow_request=req)

    def test_visualization_with_input_and_output_params_dispatches(self, tmp_path):
        """Exercise param processing: input file mapping + output filename + regular param."""
        user = _user()
        db = MagicMock()
        user_workflow = MagicMock()
        user_workflow.workflow_info = _algorithm_workflow_info()

        req = _viz_request()
        req.selectedVisualizationParams = [
            {"type": "inputFile", "name": "in", "defaultValue": "input.txt", "selectedFile": "expr.csv"},
            {"type": "outputFile", "name": "out", "defaultValue": "plot.json"},
            {"type": "string", "name": "color", "defaultValue": "red"},
        ]

        plugin_info = MagicMock()
        metadata = {"rules": {"r1": {"name": "Test Visualization"}}}

        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.get_plugin_path", return_value=("/plug/local/TestPlugin", "local")), \
             patch("app.workflow.service.validate_file_paths", return_value=["/abs/expr.csv"]), \
             patch("app.workflow.service.generate_cache_key", return_value="ck-p"), \
             patch("app.workflow.service.maybe_cleanup_cache"), \
             patch("app.workflow.service.check_cache_with_expiry", return_value=(False, None)), \
             patch("app.workflow.service.os.makedirs"), \
             patch("app.workflow.service.crud_plugin.get_plugin_by_name", return_value=plugin_info), \
             patch("builtins.open", new_callable=MagicMock), \
             patch("app.workflow.service.json.load", return_value=metadata), \
             patch("app.workflow.service.extract_rule_block", return_value=("content", "/tmp/viz/Snakefile")), \
             patch("app.workflow.service.process_data_task.apply_async", return_value=MagicMock(id="viz-p")), \
             patch("app.workflow.service.get_task_info", return_value={"task_status": "PENDING"}):
            result = service.visualize_data(db=db, user=user, workflow_request=req)

        # result filename is built from normalized params + input basename + output name
        assert result["task_id"] == "viz-p"
        assert result["result_path"] == "color_red_expr.csv_plot.json"
        assert result["cached"] is False

    def test_visualization_file_not_found_raises_domain_error(self):
        """A missing input file surfaces as the domain FileNotFoundError."""
        user = _user()
        db = MagicMock()
        user_workflow = MagicMock()
        user_workflow.workflow_info = _algorithm_workflow_info()

        req = _viz_request()
        req.selectedVisualizationParams = [
            {"type": "inputFile", "name": "in", "defaultValue": "input.txt", "selectedFile": "missing.csv"},
        ]

        from app.core.exceptions import FileNotFoundError as CellCraftFileNotFoundError
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
             patch("app.workflow.service.get_plugin_path", return_value=("/plug/local/TestPlugin", "local")), \
             patch("app.workflow.service.validate_file_paths", side_effect=FileNotFoundError("gone")):
            with pytest.raises(CellCraftFileNotFoundError):
                service.visualize_data(db=db, user=user, workflow_request=req)


# ---------------------------------------------------------------------------
# CRUD + node-data + result service functions
# ---------------------------------------------------------------------------

class TestWorkflowCrud:
    """Thin CRUD-delegating service functions."""

    def test_save_workflow_updates_existing(self):
        user = _user()
        db = MagicMock()
        wf = WorkflowCreate(id=3, title="t", thumbnail=None, workflow_info={})
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()), \
             patch("app.workflow.service.crud_workflow.update_workflow", return_value="updated") as mock_update, \
             patch("app.workflow.service.crud_workflow.create_workflow") as mock_create:
            out = service.save_workflow(db=db, user=user, workflow=wf)
        assert out == "updated"
        mock_update.assert_called_once()
        mock_create.assert_not_called()

    def test_save_workflow_creates_when_absent(self):
        user = _user()
        db = MagicMock()
        wf = WorkflowCreate(id=None, title="t", thumbnail=None, workflow_info={})
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=None), \
             patch("app.workflow.service.crud_workflow.create_workflow", return_value="created") as mock_create:
            out = service.save_workflow(db=db, user=user, workflow=wf)
        assert out == "created"
        mock_create.assert_called_once()

    def test_find_workflow_not_owned_raises_400(self):
        user = _user()
        db = MagicMock()
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=None):
            with pytest.raises(HTTPException) as exc:
                service.find_workflow(db=db, user=user, workflow_id=42)
        assert exc.value.status_code == 400
        assert "not exists" in exc.value.detail

    def test_list_workflows_empty_raises_400(self):
        user = _user()
        db = MagicMock()
        with patch("app.workflow.service.crud_workflow.get_user_workflows", return_value=[]):
            with pytest.raises(HTTPException) as exc:
                service.list_workflows(db=db, user=user)
        assert exc.value.status_code == 400

    def test_delete_workflow_not_owned_raises_400(self):
        user = _user()
        db = MagicMock()
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=None):
            with pytest.raises(HTTPException) as exc:
                service.delete_workflow(db=db, user=user, workflow_id=1)
        assert exc.value.status_code == 400

    def test_delete_workflow_removes_folder_and_db(self, tmp_path, monkeypatch):
        user = _user(username="delu")
        db = MagicMock()
        # Build a real folder inside a fake cwd so realpath guard passes.
        monkeypatch.chdir(tmp_path)
        wf_dir = tmp_path / "user" / "delu" / "workflow_1"
        wf_dir.mkdir(parents=True)
        (wf_dir / "f.txt").write_text("x")
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()), \
             patch("app.workflow.service.crud_workflow.delete_user_workflow", return_value="deleted") as mock_del:
            out = service.delete_workflow(db=db, user=user, workflow_id=1)
        assert out == "deleted"
        assert not wf_dir.exists()
        mock_del.assert_called_once()


class TestNodeData:
    """Node modal-data save/read/delete against a real tmp workflow folder."""

    def test_save_and_read_json_node(self, tmp_path, monkeypatch):
        user = _user(username="nodeu")
        db = MagicMock()
        monkeypatch.chdir(tmp_path)
        (tmp_path / "user" / "nodeu" / "workflow_1").mkdir(parents=True)

        create = WorkflowNodeFileCreate(
            id=1, node_id="n1", node_name="cfg",
            file_content=[{"k": "v"}], file_extension="json",
        )
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()):
            saved = service.save_node_data(db=db, user=user, node_file_info=create)
        assert os.path.exists(saved["file_path"])

        read = WorkflowNodeFileRead(id=1, node_id="n1", node_name="cfg", file_extension="json")
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()):
            got = service.read_node_data(db=db, user=user, node_file_info=read)
        assert got["file_content"] == [{"k": "v"}]

    def test_delete_node_removes_file(self, tmp_path, monkeypatch):
        user = _user(username="nodeu2")
        db = MagicMock()
        monkeypatch.chdir(tmp_path)
        wf_dir = tmp_path / "user" / "nodeu2" / "workflow_1"
        wf_dir.mkdir(parents=True)
        target = wf_dir / "d_n2.txt"
        target.write_text("bye")

        req = WorkflowNodeFileDelete(id=1, node_id="n2", node_name="d", file_extension="txt")
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()):
            out = service.delete_node_data(db=db, user=user, node_file_info=req)
        assert out["node_id"] == "n2"
        assert not target.exists()

    def test_read_node_missing_workflow_raises_400(self):
        user = _user()
        db = MagicMock()
        read = WorkflowNodeFileRead(id=9, node_id="n", node_name="x", file_extension="json")
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=None):
            with pytest.raises(HTTPException) as exc:
                service.read_node_data(db=db, user=user, node_file_info=read)
        assert exc.value.status_code == 400


class TestResultAccess:
    """Result listing / file access service functions."""

    def test_get_results_missing_dir_raises_404(self, tmp_path, monkeypatch):
        user = _user(username="resu")
        monkeypatch.chdir(tmp_path)
        wr = WorkflowResult(id=1, algorithm_id=1, filename=None)
        with pytest.raises(HTTPException) as exc:
            service.get_results(user=user, workflow_result=wr)
        assert exc.value.status_code == 404

    def test_get_results_lists_files(self, tmp_path, monkeypatch):
        user = _user(username="resu2")
        monkeypatch.chdir(tmp_path)
        rdir = tmp_path / "user" / "resu2" / "workflow_1" / "algorithm_1" / "results"
        rdir.mkdir(parents=True)
        (rdir / "a.csv").write_text("x,y")
        wr = WorkflowResult(id=1, algorithm_id=1, filename=None)
        out = service.get_results(user=user, workflow_result=wr)
        assert len(out) == 1
        assert out[0]["name"] == "a.csv"
        assert set(out[0].keys()) == {"name", "size", "modified_time"}

    def test_get_visualization_result_bad_file_raises_400(self):
        wr = WorkflowResult(id=1, algorithm_id=1, filename="/does/not/exist.json")
        with pytest.raises(HTTPException) as exc:
            service.get_visualization_result(workflow_result=wr)
        assert exc.value.status_code == 400

    def test_get_visualization_result_reads_json(self, tmp_path):
        f = tmp_path / "viz.json"
        f.write_text(json.dumps({"data": [], "layout": {"title": "T"}}))
        wr = WorkflowResult(id=1, algorithm_id=1, filename=str(f))
        resp = service.get_visualization_result(workflow_result=wr)
        assert isinstance(resp, JSONResponse)

    def test_check_result_returns_file_response(self, tmp_path, monkeypatch):
        user = _user(username="cru")
        monkeypatch.chdir(tmp_path)
        rdir = tmp_path / "user" / "cru" / "workflow_1" / "algorithm_1" / "results"
        rdir.mkdir(parents=True)
        (rdir / "out.csv").write_text("a,b")
        wr = WorkflowResult(id=1, algorithm_id=1, filename="out.csv")
        resp = service.check_result(user=user, workflow_result=wr)
        # The function builds a relative ./user/... path; verify it points at our file.
        assert resp.path == "./user/cru/workflow_1/algorithm_1/results/out.csv"
        assert os.path.isfile(resp.path)

    def test_check_visualization_result_reads_from_visualization_dir(self, tmp_path, monkeypatch):
        user = _user(username="cvru")
        monkeypatch.chdir(tmp_path)
        vdir = tmp_path / "user" / "cvru" / "workflow_1" / "visualization_1" / "results"
        vdir.mkdir(parents=True)
        (vdir / "viz.json").write_text(json.dumps({"data": [], "layout": {"title": "V"}}))
        wr = WorkflowResult(id=1, algorithm_id=1, filename="viz.json")
        resp = service.check_visualization_result(user=user, workflow_result=wr)
        assert isinstance(resp, JSONResponse)

    def test_check_visualization_result_missing_file_raises_400(self, tmp_path, monkeypatch):
        user = _user(username="cvru2")
        monkeypatch.chdir(tmp_path)
        vdir = tmp_path / "user" / "cvru2" / "workflow_1" / "visualization_1" / "results"
        vdir.mkdir(parents=True)  # dir exists but file does not
        wr = WorkflowResult(id=1, algorithm_id=1, filename="absent.json")
        # Avoid real sleeps during the retry loop.
        with patch("app.workflow.service.time.sleep"):
            with pytest.raises(HTTPException) as exc:
                service.check_visualization_result(user=user, workflow_result=wr)
        assert exc.value.status_code == 400
