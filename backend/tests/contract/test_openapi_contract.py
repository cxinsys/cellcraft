"""
OpenAPI contract (snapshot) test.

Freezes the API contract (paths, methods, request/response schemas) so that any
refactoring PR that unintentionally changes the public API surface fails loudly.

The snapshot itself (``openapi_snapshot.json``) is generated inside the test
environment (which runs migrations on import) via::

    python scripts/update_openapi_snapshot.py

If the contract *intentionally* changes, regenerate the snapshot with that same
command and describe the change in the PR body.

Notes:
- ``app.main`` runs ``run_migrations()`` at import time (main.py:84), so this
  test relies on the conftest ``client`` fixture, which sets ``TESTING=1`` and
  imports the app inside the configured test DB environment.
- The diff summary intentionally avoids dumping the full schema; it lists only
  the differing top-level keys and added/removed/changed paths.
"""
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

SNAPSHOT = Path(__file__).parent / "openapi_snapshot.json"

_MISSING_SNAPSHOT_MESSAGE = (
    "OpenAPI snapshot not found at {path}.\n"
    "Generate it inside the test environment (which applies migrations on import):\n"
    "    cd backend && python scripts/update_openapi_snapshot.py\n"
    "Then commit tests/contract/openapi_snapshot.json."
).format(path=SNAPSHOT)


def _summarize_diff(current: dict, expected: dict) -> str:
    """Build a concise, human-readable summary of how two OpenAPI docs differ.

    Deliberately does NOT dump the full schemas. Reports:
    - top-level keys that were added / removed / changed
    - path entries that were added / removed / changed (method-level)
    - component schema names that were added / removed / changed
    """
    lines = []

    cur_keys = set(current.keys())
    exp_keys = set(expected.keys())
    added_keys = sorted(cur_keys - exp_keys)
    removed_keys = sorted(exp_keys - cur_keys)
    changed_keys = sorted(
        k for k in (cur_keys & exp_keys) if current[k] != expected[k]
    )
    if added_keys:
        lines.append(f"  top-level keys added: {added_keys}")
    if removed_keys:
        lines.append(f"  top-level keys removed: {removed_keys}")
    if changed_keys:
        lines.append(f"  top-level keys changed: {changed_keys}")

    cur_paths = current.get("paths", {})
    exp_paths = expected.get("paths", {})
    added_paths = sorted(set(cur_paths) - set(exp_paths))
    removed_paths = sorted(set(exp_paths) - set(cur_paths))
    changed_paths = sorted(
        p for p in (set(cur_paths) & set(exp_paths)) if cur_paths[p] != exp_paths[p]
    )
    if added_paths:
        lines.append(f"  paths added ({len(added_paths)}): {added_paths}")
    if removed_paths:
        lines.append(f"  paths removed ({len(removed_paths)}): {removed_paths}")
    if changed_paths:
        lines.append(f"  paths changed ({len(changed_paths)}): {changed_paths}")

    cur_schemas = current.get("components", {}).get("schemas", {})
    exp_schemas = expected.get("components", {}).get("schemas", {})
    added_schemas = sorted(set(cur_schemas) - set(exp_schemas))
    removed_schemas = sorted(set(exp_schemas) - set(cur_schemas))
    changed_schemas = sorted(
        s for s in (set(cur_schemas) & set(exp_schemas))
        if cur_schemas[s] != exp_schemas[s]
    )
    if added_schemas:
        lines.append(f"  component schemas added ({len(added_schemas)}): {added_schemas}")
    if removed_schemas:
        lines.append(f"  component schemas removed ({len(removed_schemas)}): {removed_schemas}")
    if changed_schemas:
        lines.append(f"  component schemas changed ({len(changed_schemas)}): {changed_schemas}")

    if not lines:
        # Difference exists but not captured by the coarse categories above.
        lines.append("  (difference outside paths/components/top-level keys)")

    return "\n".join(lines)


@pytest.mark.contract
def test_openapi_schema_unchanged(client: TestClient):
    """The live OpenAPI document must match the committed snapshot exactly."""
    if not SNAPSHOT.exists():
        pytest.fail(_MISSING_SNAPSHOT_MESSAGE)

    current = client.app.openapi()
    expected = json.loads(SNAPSHOT.read_text())

    if current != expected:
        summary = _summarize_diff(current, expected)
        pytest.fail(
            "OpenAPI contract changed. If this change is intentional, regenerate "
            "the snapshot and explain the change in the PR body:\n"
            "    cd backend && python scripts/update_openapi_snapshot.py\n"
            "Differences detected:\n"
            f"{summary}"
        )
