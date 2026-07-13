"""Snapshot characterization tests for the Snakefile DAG parser (PR-10).

Purpose: pin the CURRENT parsing output of ``SnakemakeDAGParser`` so the PR-10
module split (dag_parser -> cache/parsing/status/status_analysis/status_progress)
provably preserves behavior. The expected values below were captured from the
pre-split parser and are frozen here.

These tests deliberately import via the public facade path
(``app.workflow.compiler.dag_parser``) — the same module patched by
``tests/integration/test_task_api.py`` — so they also guard the re-export.
"""
import os
import tempfile

import pytest

# Import via the facade to guard the re-export path used by callers/tests.
from app.workflow.compiler.dag_parser import SnakemakeDAGParser


REPRESENTATIVE_SNAKEFILE = '''rule TENET_Input:
    input:
        expr="user/testuser/workflow_1/algorithm_16/data/expression.csv"
    output:
        prepared="user/testuser/workflow_1/algorithm_16/results/TENET_Input.csv"
    params:
        threshold=0.05
    log:
        stdout="user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stdout",
        stderr="user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stderr"
    shell:
        "python prepare.py"

rule Run:
    input:
        prepared="user/testuser/workflow_1/algorithm_16/results/TENET_Input.csv"
    output:
        matrix="user/testuser/workflow_1/algorithm_16/results/TE_result_matrix.txt"
    log:
        stdout="user/testuser/workflow_1/algorithm_16/logs/Run.stdout",
        stderr="user/testuser/workflow_1/algorithm_16/logs/Run.stderr"
    script:
        "run_tenet.py"

rule GRN_Reconstruction__FDR:
    input:
        matrix="user/testuser/workflow_1/algorithm_16/results/TE_result_matrix.txt"
    output:
        network="user/testuser/workflow_1/algorithm_16/results/GRN.sif"
    shell:
        "python grn.py"
'''


@pytest.fixture()
def parsed_result():
    with tempfile.NamedTemporaryFile("w", suffix=".smk", delete=False) as f:
        f.write(REPRESENTATIVE_SNAKEFILE)
        path = f.name
    try:
        parser = SnakemakeDAGParser()
        yield parser.parse_snakefile_with_logs(path, use_cache=False)
    finally:
        os.unlink(path)


def test_node_ids_and_order(parsed_result):
    ids = [node["id"] for node in parsed_result["nodes"]]
    assert ids == ["TENET_Input", "Run", "GRN_Reconstruction__FDR"]


def test_node_derived_fields_snapshot(parsed_result):
    # Pin label / shortLabel / type / description derivation (frozen from pre-split parser).
    nodes = {node["id"]: node for node in parsed_result["nodes"]}

    tenet = nodes["TENET_Input"]
    assert tenet["label"] == "TENET Input"
    assert tenet["shortLabel"] == "TENET_Input"
    assert tenet["type"] == "input_processing"
    assert tenet["description"] == "Prepare input data for TENET algorithm"
    assert tenet["params"] == ["threshold"]
    assert tenet["script"] == "python prepare.py"

    run = nodes["Run"]
    assert run["label"] == "Run"
    assert run["shortLabel"] == "Run"
    assert run["type"] == "process"
    assert run["description"] == "Execute main TENET analysis"
    assert run["params"] == []
    assert run["script"] == "run_tenet.py"

    grn = nodes["GRN_Reconstruction__FDR"]
    assert grn["label"] == "GRN Reconstruction  FDR"
    assert grn["shortLabel"] == "GRN_Reconstruction\n_FDR"
    assert grn["type"] == "network_analysis"
    assert grn["description"] == "Reconstruct gene regulatory network using FDR threshold"
    assert grn["script"] == "python grn.py"


def test_node_io_and_logs_snapshot(parsed_result):
    nodes = {node["id"]: node for node in parsed_result["nodes"]}

    assert nodes["TENET_Input"]["inputs"] == [
        "user/testuser/workflow_1/algorithm_16/data/expression.csv"
    ]
    assert nodes["TENET_Input"]["outputs"] == [
        "user/testuser/workflow_1/algorithm_16/results/TENET_Input.csv"
    ]
    assert nodes["TENET_Input"]["log_paths"] == {
        "stdout": "user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stdout",
        "stderr": "user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stderr",
    }
    # No ``log:`` section -> empty log_paths (frozen behavior).
    assert nodes["GRN_Reconstruction__FDR"]["log_paths"] == {}


def test_edges_snapshot(parsed_result):
    assert parsed_result["edges"] == [
        {"source": "TENET_Input", "target": "Run", "label": "TENET_Input.csv"},
        {"source": "Run", "target": "GRN_Reconstruction__FDR", "label": "TE_result_matrix.txt"},
    ]


def test_execution_sequence_snapshot(parsed_result):
    assert parsed_result["execution_sequence"] == [
        "TENET_Input",
        "Run",
        "GRN_Reconstruction__FDR",
    ]


def test_empty_snakefile_returns_empty_graph():
    with tempfile.NamedTemporaryFile("w", suffix=".smk", delete=False) as f:
        f.write("# just a comment, no rules\n")
        path = f.name
    try:
        result = SnakemakeDAGParser().parse_snakefile_with_logs(path, use_cache=False)
        assert result == {"nodes": [], "edges": [], "execution_sequence": []}
    finally:
        os.unlink(path)
