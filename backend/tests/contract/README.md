# Contract & characterization test guide

This directory holds the **API contract** safety net for the v1.0.7 refactoring
(PR-1, phase 0). It exists to mechanically prove that refactoring PRs do **not**
change the public API surface or observable behavior.

## What lives here

- `test_openapi_contract.py` — compares the live `app.openapi()` document against
  the committed snapshot `openapi_snapshot.json`. Any change to paths, methods,
  or request/response schemas fails the test with a concise diff summary.
- `openapi_snapshot.json` — the frozen OpenAPI document. **Generated inside the
  test environment**, not by hand (see below).

Characterization tests that pin *behavior* (not just the schema) live alongside
the existing integration tests in `backend/tests/integration/`:

- `test_characterization_auth.py`
- `test_characterization_plugin_upload.py`
- `test_characterization_workflow_compile.py`
- `test_characterization_task_status.py`
- `test_characterization_file_upload.py`

## Generating / updating the OpenAPI snapshot

`app.main` runs Alembic migrations at import time (main.py:84), so the snapshot
**must** be generated where a database is reachable. Use the test compose stack:

```bash
# 1. Build the backend app image once (arm64 example; pick your platform/Dockerfile),
#    then the test image (app image + test deps from requirements-test.txt).
docker build -f backend/Dockerfile.cpu \
  --build-arg TARGETPLATFORM=linux/arm64 --platform linux/arm64 \
  -t cellcraft-backend:refactor-test ./backend
docker compose -f docker-compose.test.yml build backend-test

# 2. Start the isolated test stack (own network + own DB).
docker compose -f docker-compose.test.yml up -d

# 3. Generate the snapshot inside the container.
docker compose -f docker-compose.test.yml exec backend-test \
    python scripts/update_openapi_snapshot.py

# 4. The snapshot is written to backend/tests/contract/openapi_snapshot.json
#    (the ./backend bind mount surfaces it on the host). Commit it.
```

Update the snapshot **only** when an API change is intentional, and describe the
change in the PR body.

## Running the tests

The test suite talks to a PostgreSQL test database. The connection string is read
from the `TEST_DATABASE_URI` environment variable (falling back to
`postgresql://test_user:test_pass@localhost:5433/cellcraft_test`).

### Using the test compose override (recommended)

```bash
# Build images (see above) + bring up test-db and backend-test.
docker compose -f docker-compose.test.yml up -d

# Run the whole suite inside the container (TESTING=1 and TEST_DATABASE_URI are
# already set in the compose file, pointing at the in-network test-db).
docker compose -f docker-compose.test.yml exec backend-test pytest

# Run just the contract + characterization tests.
docker compose -f docker-compose.test.yml exec backend-test \
    pytest tests/contract tests/integration -m "contract or characterization"

# Tear down when done.
docker compose -f docker-compose.test.yml down
```

If host port 5433 is occupied, override it:

```bash
TEST_DB_PORT=5544 docker compose -f docker-compose.test.yml up -d
```

(The `backend-test` service reaches the DB over the internal network on port
5432, so overriding the host port does not affect in-container test runs.)

### Running against a local/dev test-db

The existing `docker-compose.dev.yml` also defines a `test-db` service on host
port 5433. To run the suite from a local checkout against it:

```bash
docker compose -f docker-compose.dev.yml up test-db -d
cd backend && pytest
```

To point at a different host/port, set `TEST_DATABASE_URI` (and `TEST_DB_PORT`
for compose) accordingly.
