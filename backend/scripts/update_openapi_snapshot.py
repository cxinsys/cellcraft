#!/usr/bin/env python
"""
Generate / refresh the OpenAPI contract snapshot.

This MUST run inside the test environment. ``app.main`` executes
``run_migrations()`` at import time (main.py:84), so importing the app requires a
reachable database. Setting ``TESTING=1`` before the import routes the app to the
test DB configuration (see ``app/db/session.py`` and conftest).

Usage (from the ``backend`` directory, inside the container / test env):

    python scripts/update_openapi_snapshot.py

The snapshot is written with ``sort_keys=True`` and ``indent=2`` so diffs stay
stable and reviewable.
"""
import json
import os
import sys
from pathlib import Path

# Ensure the test DB configuration is selected BEFORE importing the app,
# mirroring the conftest import pattern (main.py runs migrations on import).
os.environ.setdefault("TESTING", "1")

# Make sure the backend package root is importable when run as a script.
BACKEND_ROOT = Path(__file__).resolve().parent.parent
if str(BACKEND_ROOT) not in sys.path:
    sys.path.insert(0, str(BACKEND_ROOT))

SNAPSHOT_PATH = BACKEND_ROOT / "tests" / "contract" / "openapi_snapshot.json"


def generate_snapshot() -> Path:
    """Import the app, dump its OpenAPI document, and write the snapshot file."""
    from app.main import app  # imported here so TESTING is already set

    openapi_schema = app.openapi()

    SNAPSHOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SNAPSHOT_PATH.write_text(
        json.dumps(openapi_schema, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    )
    return SNAPSHOT_PATH


if __name__ == "__main__":
    written = generate_snapshot()
    print(f"OpenAPI snapshot written to: {written}")
