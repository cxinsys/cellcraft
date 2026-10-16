"""
Golden Snakefile tests for official plugins (#137).

``tests/fixtures/snakefile_golden/expected`` was generated with the pre-#137
implementation (release/v1.0.8 a57c477) by running frontend-shaped requests
through ``compile_workflow`` / ``visualize_data`` (see ``build_cases.py`` and
``generate_golden.py``). Input validation, placeholder filtering and literal
escaping must not change the generated Snakefile for these requests.
"""
import os
import sys

import pytest

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures", "snakefile_golden")
sys.path.insert(0, GOLDEN_DIR)
import harness  # noqa: E402

CASES = harness.load_cases()


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_service_generated_snakefile_matches_golden(case, tmp_path):
    generated = harness.render(case, str(tmp_path))

    with open(os.path.join(GOLDEN_DIR, "expected", f"{case['id']}.Snakefile"), encoding="utf-8") as f:
        assert generated == f.read()


def test_all_official_plugins_covered():
    plugins = {case["plugin"] for case in CASES}
    assert plugins == {"FastSCODE", "FastTENET", "GENIE3", "GRNBoost2", "GRNViz", "LEAP", "Scribe", "TENET"}
    assert {case["id"] for case in CASES if case["kind"] == "visualization"} == {
        "GRNViz_Heatmap", "GRNViz_Barplot", "GRNViz_Network"}
