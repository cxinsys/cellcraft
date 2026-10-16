"""
Golden Snakefile tests for official plugins (#137).

``tests/fixtures/snakefile_golden/expected`` was generated with the pre-#137
implementation (release/v1.0.8 a57c477) from the official plugin templates
(release/plugins-v1.1) and the values the frontend sends. Input validation and
literal escaping must not change the generated Snakefile for these values.
"""
import json
import os

import pytest

from app.workflow.compiler.snakefile import change_snakefile_parameter
from app.workflow.compiler.visualization import extract_rule_block, format_visualization_param

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures", "snakefile_golden")
VISUALIZATION_CONTEXT_KEYS = {
    "user_name", "workflow_id", "algorithm_id", "visualization_id",
    "plugin_name", "visualization_name", "result_filename",
}

with open(os.path.join(GOLDEN_DIR, "cases.json"), encoding="utf-8") as _f:
    CASES = json.load(_f)


def _case_id(case):
    if case["kind"] == "visualization":
        return f"{case['plugin']}_{case['rule'].replace(' ', '_')}"
    return case["plugin"]


@pytest.mark.parametrize("case", CASES, ids=_case_id)
def test_generated_snakefile_matches_golden(case, tmp_path):
    template = os.path.join(GOLDEN_DIR, "templates", case["plugin"], "Snakefile")
    output = str(tmp_path / "Snakefile")

    if case["kind"] == "analysis":
        change_snakefile_parameter(template, output, dict(case["user_input"]), dict(case["plugin_params"]))
    else:
        # visualize_data와 동일하게 사용자 파라미터만 리터럴 변환 후 시스템 placeholder와 병합
        params = {
            name: (value if name in VISUALIZATION_CONTEXT_KEYS else format_visualization_param(value))
            for name, value in case["parameters"].items()
        }
        extract_rule_block(
            snakefile_path=template, visualization_title=case["rule"], result_path=output,
            parameters=params, file_mappings=dict(case["file_mappings"]),
        )

    with open(output, encoding="utf-8") as f:
        generated = f.read()
    with open(os.path.join(GOLDEN_DIR, "expected", f"{_case_id(case)}.Snakefile"), encoding="utf-8") as f:
        expected = f.read()
    assert generated == expected


def test_all_official_plugins_covered():
    plugins = {case["plugin"] for case in CASES}
    assert plugins == {"FastSCODE", "FastTENET", "GENIE3", "GRNBoost2", "GRNViz", "LEAP", "Scribe", "TENET"}
