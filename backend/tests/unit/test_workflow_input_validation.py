"""
Server-side validation of client-supplied workflow values (#137).

Covers the validators, Snakefile literal formatting, and the service entry
points that build filesystem paths from frontend input (compile, visualization,
result download, node data). Legitimate frontend values (see module docstring of
``app.workflow.validators``) must keep working; malicious values must be
rejected before any filesystem change.
"""
import ast
import json
import os
from unittest.mock import MagicMock, patch

import pytest
from fastapi import HTTPException

from app.core.exceptions import (
    NotFoundError, PermissionDeniedError, ValidationError, ValidationFailedError,
)
from app.workflow import service
from app.workflow import validators as v
from app.workflow.compiler.snakefile import change_snakefile_parameter, format_plugin_param
from app.workflow.compiler.visualization import format_visualization_param
from app.workflow.schemas import (
    WorkflowCreate, WorkflowNodeFileCreate, WorkflowNodeFileRead, WorkflowNodeFileDelete,
    WorkflowResult, WorkflowVisualizationRequest,
)


def _user(username="john", user_id=7):
    u = MagicMock()
    u.id = user_id
    u.username = username
    return u


# ---------------------------------------------------------------------------
# validators
# ---------------------------------------------------------------------------

class TestNodeId:
    @pytest.mark.parametrize("value", [3, "3", "12", "node_1", "a-b"])
    def test_frontend_ids_accepted(self, value):
        assert v.validate_node_id(value) == str(value)

    @pytest.mark.parametrize("value", ["../2", "2/../../x", "", None, True, "a" * 65, "1 2", "1.2"])
    def test_rejected(self, value):
        with pytest.raises(ValidationFailedError):
            v.validate_node_id(value)


class TestBareFilename:
    @pytest.mark.parametrize("value", [
        "pbmc_3k.h5ad",                     # sanitized upload
        "pbmc light 1000.h5ad",             # shared tutorial file (on-disk name)
        "인간_TF.txt",
        "Top_Genes_10_Sample_Size_100_expression.csv_Heatmap.json",  # visualization result
        "GRN.sif",
    ])
    def test_frontend_names_accepted(self, value):
        assert v.validate_bare_filename(value) == value

    @pytest.mark.parametrize("value", [
        "../secret.txt", "a/b.csv", "..", ".", "", None, 3,
        "a\\b.csv", 'x".h5ad', "x\n.h5ad", "x\x00.h5ad", "a" * 256,
    ])
    def test_rejected(self, value):
        with pytest.raises(ValidationFailedError):
            v.validate_bare_filename(value)


class TestNodeData:
    def test_scatter_plot_values_accepted(self):
        assert v.validate_node_name("ScatterPlot") == "ScatterPlot"
        assert v.validate_node_data_extension("json") == "json"

    @pytest.mark.parametrize("value", ["../x", "a/b", "Scatter Plot", ""])
    def test_bad_node_name(self, value):
        with pytest.raises(ValidationFailedError):
            v.validate_node_name(value)

    @pytest.mark.parametrize("value", ["py", "json/../x", "", None])
    def test_bad_extension(self, value):
        with pytest.raises(ValidationFailedError):
            v.validate_node_data_extension(value)


class TestRelativeOutputPath:
    @pytest.mark.parametrize("value, expected", [
        ("GRN.sif", "GRN.sif"), ("plots/out.png", "plots/out.png"), ("./out.csv", "out.csv"),
        ("Top_Genes_10_expression.csv_Heatmap.json", "Top_Genes_10_expression.csv_Heatmap.json"),
    ])
    def test_accepted(self, value, expected):
        assert v.validate_relative_output_path(value) == expected

    @pytest.mark.parametrize("value", ["../x", "a/../../x", "/abs/x", "..", ".", "", None, 'a"b', "a\\b"])
    def test_rejected(self, value):
        with pytest.raises(ValidationFailedError):
            v.validate_relative_output_path(value)

    def test_sanitize_filename_part(self):
        assert v.sanitize_filename_part('Group_a/b"c') == "Group_a_b_c"
        assert v.sanitize_filename_part("Top_Genes_10") == "Top_Genes_10"


class TestPluginPlaceholders:
    def _write_metadata(self, tmp_path, rules):
        (tmp_path / "metadata.json").write_text(json.dumps({"rules": rules}))

    RULES = {
        "0": {"name": "Input", "parameters": [
            {"name": "input", "type": "inputFile", "defaultValue": "input.h5ad"},
            {"name": "geneList", "type": "optionalInputFile", "defaultValue": "geneList.txt"},
            {"name": "cell group", "type": "h5adParameter", "defaultValue": ""},
            {"name": "expression", "type": "outputFile", "defaultValue": "expression.csv"},
        ]},
        "1": {"name": "Run", "parameters": [
            {"name": "number of devices", "type": "int", "defaultValue": "1"},
            {"name": "log", "type": "boolean", "defaultValue": "false"},
        ]},
    }

    def test_all_rules(self, tmp_path):
        self._write_metadata(tmp_path, self.RULES)
        assert v.load_plugin_placeholders(str(tmp_path)) == {
            "params": {"cell group", "number of devices", "log"},
            "input_files": {"input.h5ad", "geneList.txt"},
        }

    def test_single_rule(self, tmp_path):
        self._write_metadata(tmp_path, self.RULES)
        assert v.load_plugin_placeholders(str(tmp_path), "Run") == {
            "params": {"number of devices", "log"}, "input_files": set()}

    def test_rules_as_list(self, tmp_path):
        self._write_metadata(tmp_path, list(self.RULES.values()))
        assert v.load_plugin_placeholders(str(tmp_path))["input_files"] == {"input.h5ad", "geneList.txt"}

    @pytest.mark.parametrize("rules", ["x", 3, {"0": "not a rule"}, {"0": {"parameters": "x"}}, [None]])
    def test_malformed_metadata_does_not_raise(self, tmp_path, rules):
        self._write_metadata(tmp_path, rules)
        result = v.load_plugin_placeholders(str(tmp_path))
        assert result is None or result == {"params": set(), "input_files": set()}

    def test_missing_metadata(self, tmp_path):
        assert v.load_plugin_placeholders(str(tmp_path)) is None
        assert v.placeholder_names(None, "params") is None

    def test_filter_by_metadata(self):
        params = {"cell group": "leiden", "log": "false", "params.FDR": "x", "unknown": 1}
        assert v.filter_plugin_params(params, {"cell group", "log"}) == {"cell group": "leiden", "log": "false"}

    def test_fallback_excludes_snakemake_namespaces(self):
        params = {"cell group": "leiden", "log": "false", "params.FDR": "x", "input.input": "y", "log.stdout": "z"}
        assert v.filter_plugin_params(params, None) == {"cell group": "leiden", "log": "false"}


# ---------------------------------------------------------------------------
# Snakefile literal formatting
# ---------------------------------------------------------------------------

class TestFormatPluginParam:
    @pytest.mark.parametrize("value, expected", [
        ("leiden", '"leiden"'),
        ("10", '"10"'),
        (10, "10"),
        (0.01, "0.01"),
        (True, "True"),
        (["0", "1", "B cell"], '["0", "1", "B cell"]'),
        ([1, 2], "[1, 2]"),
        ({"a": "x", "b": 1}, '{a: "x", b: 1}'),
        ("한글 라벨", '"한글 라벨"'),
    ])
    def test_unchanged_for_normal_values(self, value, expected):
        assert format_plugin_param(value) == expected

    @pytest.mark.parametrize("value", [
        'x", __import__("os").system("id"), "',
        "a\\", "line\nbreak",
        ['ok', '"]; import os #'],
    ])
    def test_escaped_value_round_trips_as_literal(self, value):
        assert ast.literal_eval(format_plugin_param(value)) == value

    def test_dict_key_must_be_identifier(self):
        with pytest.raises(ValueError):
            format_plugin_param({"a: 1, b": "x"})

    def test_unsupported_type(self):
        with pytest.raises(ValueError):
            format_plugin_param(None)

    def test_malicious_param_cannot_add_statements(self, tmp_path):
        template = tmp_path / "Snakefile"
        template.write_text("rule a:\n    params:\n        cell_group={cell group},\n")
        out = tmp_path / "out"
        payload = 'x", os.system("id"), "'
        change_snakefile_parameter(str(template), str(out), {}, {"cell group": payload})

        line = out.read_text().splitlines()[2].strip()
        literal = line[len("cell_group="):-1]
        assert ast.literal_eval(literal) == payload


class TestFormatVisualizationParam:
    @pytest.mark.parametrize("value, expected", [
        ("10", "10"), (" 20 ", " 20 "), ("-1.5", "-1.5"), ("1e-3", "1e-3"),
        (10, "10"), (0.5, "0.5"), (True, "True"), (None, "None"),
        (["a", 1], "['a', 1]"),
    ])
    def test_previous_str_output_kept(self, value, expected):
        assert format_visualization_param(value) == expected

    @pytest.mark.parametrize("value", ["leiden", '1); import os; os.system("id"', "a\nb"])
    def test_non_numeric_strings_become_literals(self, value):
        assert ast.literal_eval(format_visualization_param(value)) == value


# ---------------------------------------------------------------------------
# compile_workflow
# ---------------------------------------------------------------------------

def _plugin_dir(tmp_path):
    plugin = tmp_path / "plugin" / "P"
    plugin.mkdir(parents=True)
    (plugin / "Snakefile").write_text(
        'rule a:\n    input:\n        input="{input.h5ad}"\n'
        '    output:\n        out="user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results/out.csv"\n'
        '    params:\n        cell_group={cell group}\n'
        '    shell:\n        "run {input.input} {params.cell_group} > {output.out}"\n'
    )
    (plugin / "metadata.json").write_text(json.dumps({"rules": {"0": {"name": "a", "parameters": [
        {"name": "input", "type": "inputFile", "defaultValue": "input.h5ad"},
        {"name": "cell group", "type": "h5adParameter", "defaultValue": ""},
        {"name": "out", "type": "outputFile", "defaultValue": "out.csv"},
    ]}}}))
    return plugin


def _algorithm_node(node_id=3, file_name="pbmc_3k.h5ad", output="out.csv", params=None,
                    input_key="input.h5ad", outputs=None):
    rules = [{"parameters": params or [{"name": "cell group", "defaultValue": "leiden"}]}]
    return {
        "id": node_id, "name": "algorithm", "class": "Algorithm", "outputs": outputs or {},
        "data": {
            "selectedPlugin": {"name": "P", "source": "local"},
            "selectedPluginInputOutput": [
                {"type": "inputFile", "defaultValue": input_key, "file_name": file_name, "fileSource": "user"},
                {"type": "outputFile", "defaultValue": output, "activate": True},
            ],
            "selectedPluginRules": rules,
        },
    }


def _workflow_info(node_id=3, **kwargs):
    return {"drawflow": {"Home": {"data": {str(node_id): _algorithm_node(node_id, **kwargs)}}}}


def _compile(tmp_path, monkeypatch, workflow_info):
    monkeypatch.chdir(tmp_path)
    plugin = _plugin_dir(tmp_path)
    user_workflow = MagicMock()
    user_workflow.workflow_info = workflow_info
    db = MagicMock()
    db.query.return_value.filter_by.return_value.first.return_value = None
    workflow = WorkflowCreate(id=1, title="t", thumbnail=None, workflow_info=workflow_info)
    dispatch = MagicMock(return_value=MagicMock(id="tid-1"))
    with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=user_workflow), \
         patch("app.workflow.service.crud_workflow.update_workflow") as update, \
         patch("app.workflow.service.get_plugin_path", return_value=(str(plugin), "local")), \
         patch("app.workflow.service.process_data_task.apply_async", dispatch), \
         patch("app.workflow.service.get_task_info", return_value={"task_status": "PENDING"}):
        try:
            result = service.compile_workflow(db=db, user=_user(), workflow=workflow)
        finally:
            _compile.update = update
    return result, dispatch


class TestCompileWorkflow:
    def test_frontend_workflow_compiles(self, tmp_path, monkeypatch):
        result, dispatch = _compile(tmp_path, monkeypatch, _workflow_info())

        assert result["task_ids"] == ["tid-1"]
        snakefile = (tmp_path / "user/john/workflow_1/algorithm_3/Snakefile").read_text()
        assert 'input="user/john/data/pbmc_3k.h5ad"' in snakefile
        assert 'cell_group="leiden"' in snakefile
        assert dispatch.call_args.kwargs["args"][3] == ["./user/john/workflow_1/algorithm_3/results/out.csv"]

    def test_undeclared_param_not_substituted(self, tmp_path, monkeypatch):
        params = [{"name": "cell group", "defaultValue": "leiden"},
                  {"name": "params.cell_group", "defaultValue": "x; rm -rf /"}]
        _compile(tmp_path, monkeypatch, _workflow_info(params=params))

        snakefile = (tmp_path / "user/john/workflow_1/algorithm_3/Snakefile").read_text()
        assert "{params.cell_group}" in snakefile
        assert "rm -rf" not in snakefile

    def test_output_subdirectory_allowed(self, tmp_path, monkeypatch):
        _, dispatch = _compile(tmp_path, monkeypatch, _workflow_info(output="plots/out.csv"))
        assert dispatch.call_args.kwargs["args"][3] == ["./user/john/workflow_1/algorithm_3/results/plots/out.csv"]

    def test_undeclared_input_placeholder_not_substituted(self, tmp_path, monkeypatch):
        _compile(tmp_path, monkeypatch, _workflow_info(input_key="output.out"))
        snakefile = (tmp_path / "user/john/workflow_1/algorithm_3/Snakefile").read_text()
        assert "{output.out}" in snakefile

    def test_invalid_second_algorithm_blocks_whole_compile(self, tmp_path, monkeypatch):
        data = {"3": _algorithm_node(3), "4": _algorithm_node(4, file_name="../../other/x.h5ad")}
        info = {"drawflow": {"Home": {"data": data}}}
        dispatch = None
        with pytest.raises(HTTPException) as exc:
            _, dispatch = _compile(tmp_path, monkeypatch, info)

        assert exc.value.status_code == 400
        assert dispatch is None
        _compile.update.assert_not_called()
        assert not (tmp_path / "user").exists()

    def test_connected_visualization_id_validated(self, tmp_path, monkeypatch):
        vis_key = "9/../../../other"
        outputs = {"output_1": {"connections": [{"node": vis_key, "output": "input_1"}]}}
        info = _workflow_info(outputs=outputs)
        info["drawflow"]["Home"]["data"][vis_key] = {"id": vis_key, "class": "Visualization", "data": {}}
        with pytest.raises(HTTPException) as exc:
            _compile(tmp_path, monkeypatch, info)
        assert exc.value.status_code == 400
        assert not (tmp_path / "user").exists()

    @pytest.mark.parametrize("info", [
        _workflow_info(node_id="../../other"),
        _workflow_info(file_name="../../other/data/x.h5ad"),
        _workflow_info(output="../../../escape.csv"),
    ], ids=["node_id", "input_file", "output_file"])
    def test_rejected_before_filesystem_change(self, tmp_path, monkeypatch, info):
        with pytest.raises(HTTPException) as exc:
            _compile(tmp_path, monkeypatch, info)

        assert exc.value.status_code == 400
        assert not (tmp_path / "user").exists()
        assert not (tmp_path / "escape.csv").exists()


# ---------------------------------------------------------------------------
# visualization
# ---------------------------------------------------------------------------

def _visualization_request(**overrides):
    payload = dict(
        id=1, current_node_id="7", algorithm_id="3",
        selectedVisualizationPlugin={"name": "GRNViz", "source": "official"},
        selectedScript={"name": "Barplot"},
        selectedVisualizationParams=[
            {"type": "inputFile", "name": "input", "defaultValue": "input.txt", "selectedFile": "GRN.sif"},
            {"type": "int", "name": "Top Genes", "defaultValue": "10"},
            {"type": "outputFile", "name": "Barplot", "defaultValue": "Barplot.json"},
        ],
    )
    payload.update(overrides)
    return WorkflowVisualizationRequest(**payload)


class TestVisualizeDataValidation:
    @pytest.mark.parametrize("overrides, field", [
        ({"current_node_id": "../../x"}, "current_node_id"),
        ({"algorithm_id": "3/../../x"}, "algorithm_id"),
        ({"selectedVisualizationParams": [
            {"type": "inputFile", "name": "input", "defaultValue": "input.txt", "selectedFile": "../../x.sif"}]},
         "selectedFile"),
        ({"selectedVisualizationParams": [
            {"type": "outputFile", "name": "Barplot", "defaultValue": "../../escape.json"}]},
         "outputFile"),
    ])
    def test_rejected(self, tmp_path, monkeypatch, overrides, field):
        monkeypatch.chdir(tmp_path)
        db = MagicMock()
        with patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock()), \
             patch("app.workflow.service.get_plugin_path", return_value=(str(tmp_path), "official")), \
             patch("app.workflow.service.validate_file_paths", return_value=[str(tmp_path / "GRN.sif")]), \
             patch("app.workflow.service.maybe_cleanup_cache"), \
             patch("app.workflow.service.check_cache_with_expiry", return_value=(False, None)):
            with pytest.raises(ValidationError) as exc:
                service.visualize_data(db=db, user=_user(), workflow_request=_visualization_request(**overrides))

        assert exc.value.status_code == 400
        assert field in json.dumps(exc.value.detail)
        assert not (tmp_path / "user").exists()


# ---------------------------------------------------------------------------
# results & node data
# ---------------------------------------------------------------------------

class TestResultFiles:
    @pytest.fixture
    def results_dir(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        results = tmp_path / "user/john/workflow_1/algorithm_3/results"
        results.mkdir(parents=True)
        (results / "GRN.sif").write_text("a b")
        (tmp_path / "user/other").mkdir(parents=True)
        (tmp_path / "user/other/secret.txt").write_text("secret")
        return results

    def test_listed_file_served(self, results_dir):
        response = service.check_result(user=_user(), workflow_result=WorkflowResult(id=1, algorithm_id=3, filename="GRN.sif"))
        assert response.path.endswith("results/GRN.sif")

    def test_traversal_filename_rejected(self, results_dir):
        with pytest.raises(ValidationFailedError):
            service.check_result(user=_user(), workflow_result=WorkflowResult(
                id=1, algorithm_id=3, filename="../../../other/secret.txt"))

    def test_missing_file_404(self, results_dir):
        with pytest.raises(NotFoundError):
            service.check_result(user=_user(), workflow_result=WorkflowResult(id=1, algorithm_id=3, filename="missing.sif"))

    def test_visualization_result_traversal_rejected(self, results_dir):
        with pytest.raises(ValidationFailedError):
            service.check_visualization_result(
                user=_user(), workflow_result=WorkflowResult(id=1, algorithm_id=3, filename="../../../other/secret.txt"))

    def test_visualize_result_restricted_to_user_dir(self, results_dir, tmp_path):
        own = tmp_path / "user/john/workflow_1/plot.json"
        own.write_text('{"data": []}')
        response = service.get_visualization_result(
            user=_user(), workflow_result=WorkflowResult(filename="./user/john/workflow_1/plot.json"))
        assert json.loads(response.body) == {"data": []}

        with pytest.raises(PermissionDeniedError):
            service.get_visualization_result(
                user=_user(), workflow_result=WorkflowResult(filename="./user/other/secret.txt"))


class TestNodeDataFiles:
    @pytest.fixture(autouse=True)
    def _cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "user/john/workflow_1").mkdir(parents=True)

    def _patch_owner(self):
        return patch("app.workflow.service.crud_workflow.get_user_workflow", return_value=MagicMock())

    def test_scatter_plot_round_trip(self, tmp_path):
        info = dict(id=1, node_id="4", node_name="ScatterPlot", file_extension="json")
        with self._patch_owner():
            saved = service.save_node_data(db=MagicMock(), user=_user(),
                                           node_file_info=WorkflowNodeFileCreate(file_content=[1, 2], **info))
            read = service.read_node_data(db=MagicMock(), user=_user(), node_file_info=WorkflowNodeFileRead(**info))
            service.delete_node_data(db=MagicMock(), user=_user(), node_file_info=WorkflowNodeFileDelete(**info))

        assert saved["file_path"] == "./user/john/workflow_1/ScatterPlot_4.json"
        assert read["file_content"] == [1, 2]
        assert not (tmp_path / "user/john/workflow_1/ScatterPlot_4.json").exists()

    @pytest.mark.parametrize("info", [
        dict(node_id="4", node_name="../../../escape", file_extension="json"),
        dict(node_id="4/../../../escape", node_name="ScatterPlot", file_extension="json"),
        dict(node_id="4", node_name="ScatterPlot", file_extension="json/../../../escape"),
    ])
    def test_traversal_rejected(self, tmp_path, info):
        with self._patch_owner(), pytest.raises(ValidationFailedError):
            service.save_node_data(db=MagicMock(), user=_user(),
                                   node_file_info=WorkflowNodeFileCreate(id=1, file_content=[1], **info))
        assert not any(p.name.startswith("escape") for p in tmp_path.rglob("*"))


# ---------------------------------------------------------------------------
# plugin container command
# ---------------------------------------------------------------------------

class TestSnakemakeCommandQuoting:
    def test_targets_and_paths_shell_quoted(self, tmp_path, monkeypatch):
        import contextlib
        import shlex
        from app.workflow.compiler import snakefile as sf

        monkeypatch.chdir(tmp_path)
        monkeypatch.setenv("HOSTNAME", "celery-container")
        (tmp_path / "user/john/workflow_1/visualization_7").mkdir(parents=True)

        client = MagicMock()
        client.containers.get.return_value.attrs = {"Mounts": [{"Destination": "/app", "Source": "/srv/backend"}]}
        client.containers.run.return_value = MagicMock(id="cid", status="running")
        client.api.exec_create.return_value = {"Id": "exec-id"}
        client.api.exec_start.return_value = iter([b"done"])
        client.api.exec_inspect.return_value = {"ExitCode": 0}
        db = MagicMock()
        db.query.return_value.filter_by.return_value.first.return_value = None

        @contextlib.contextmanager
        def fake_session():
            yield db

        target = "./user/john/workflow_1/visualization_7/results/Top Genes_10; touch pwned_Barplot.json"
        with patch.object(sf.docker, "from_env", return_value=client), \
             patch("app.db.session.get_db_session", fake_session), \
             patch.object(sf, "wait_for_container_ready"), \
             patch.object(sf.container_manager, "register_container"), \
             patch.object(sf.container_manager, "unregister_container"), \
             patch.object(sf.container_manager, "_is_cleanup_in_progress", return_value=False):
            sf.exec_in_plugin("P", "./user/john/workflow_1/visualization_7/Snakefile", [target], str(tmp_path))

        bash_cmd = client.api.exec_create.call_args.args[1][2]
        assert shlex.quote(target) in bash_cmd
        tokens = shlex.split(bash_cmd)
        assert target in tokens
        assert "touch" not in tokens


# ---------------------------------------------------------------------------
# visualization through the service (golden harness)
# ---------------------------------------------------------------------------

def _barplot_case(params):
    import copy
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures", "snakefile_golden"))
    import harness
    case = copy.deepcopy(next(c for c in harness.load_cases() if c["id"] == "GRNViz_Barplot"))
    case["request"]["selectedVisualizationParams"] = params
    return harness, case


class TestVisualizationResultName:
    INPUT = {"type": "inputFile", "name": "input", "defaultValue": "input.txt", "selectedFile": "Numlinks.sif"}

    def test_param_value_with_separator_sanitized(self, tmp_path):
        harness, case = _barplot_case([
            self.INPUT,
            {"type": "int", "name": "Top Genes", "defaultValue": "10/../x"},
            {"type": "outputFile", "name": "Barplot", "defaultValue": "Barplot.json"},
        ])
        snakefile = harness.render(case, str(tmp_path))
        assert 'results/Top_Genes_10_.._x_Numlinks.sif_Barplot.json"' in snakefile

    def test_empty_output_uses_default_name(self, tmp_path):
        harness, case = _barplot_case([
            self.INPUT,
            {"type": "int", "name": "Top Genes", "defaultValue": "10"},
            {"type": "outputFile", "name": "Barplot", "defaultValue": ""},
        ])
        snakefile = harness.render(case, str(tmp_path))
        assert 'results/Top_Genes_10_Numlinks.sif_visualization_result.json"' in snakefile

    def test_undeclared_file_mapping_ignored(self, tmp_path):
        harness, case = _barplot_case([
            self.INPUT,
            {"type": "inputFile", "name": "x", "defaultValue": "log.stdout", "selectedFile": "Numlinks.sif"},
            {"type": "int", "name": "Top Genes", "defaultValue": "10"},
            {"type": "outputFile", "name": "Barplot", "defaultValue": "Barplot.json"},
        ])
        snakefile = harness.render(case, str(tmp_path))
        assert "{log.stdout}" in snakefile
