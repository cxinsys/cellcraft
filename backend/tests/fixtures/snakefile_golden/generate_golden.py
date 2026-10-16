"""Regenerate ``expected/`` from a reference revision (see build_cases.py)."""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CELLCRAFT_BACKEND", os.path.abspath(os.path.join(HERE, "..", "..", ".."))))
sys.path.insert(0, HERE)
for key, value in {"TESTING": "1", "POSTGRES_USER": "x", "POSTGRES_PASSWORD": "x",
                   "POSTGRES_HOST": "x", "POSTGRES_PORT": "1", "POSTGRES_DB": "x"}.items():
    os.environ.setdefault(key, value)

import harness  # noqa: E402

os.makedirs(os.path.join(HERE, "expected"), exist_ok=True)
for case in harness.load_cases():
    with tempfile.TemporaryDirectory() as workdir:
        content = harness.render(case, workdir)
    with open(os.path.join(HERE, "expected", f"{case['id']}.Snakefile"), "w", encoding="utf-8") as f:
        f.write(content)
    print(case["id"])
