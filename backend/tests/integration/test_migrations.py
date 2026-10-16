"""
Startup migration scenarios (#138: first revision after the baseline).

Each scenario runs ``run_migrations`` against its own temporary database on the
test PostgreSQL server:

- A: empty database (new install)
- B: database created by ``create_all()`` before Alembic (pre-v1.0.7 install)
- C: Alembic-managed database at the baseline revision (v1.0.7 install)
- D: version table left without application tables
"""
import os
import uuid

import pytest
from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.engine import make_url

from app.core.startup import BASELINE_REVISION, run_migrations

HEAD_REVISION = "0002_must_change_password"
BACKEND_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


@pytest.fixture
def temp_database_url(monkeypatch):
    monkeypatch.chdir(BACKEND_ROOT)  # alembic.ini / alembic/ are resolved relative to the backend root
    base_url = make_url(os.environ.get(
        "TEST_DATABASE_URI", "postgresql://test_user:test_pass@localhost:5433/cellcraft_test"))
    name = f"cellcraft_migration_{uuid.uuid4().hex[:8]}"
    admin_engine = create_engine(base_url.set(database="postgres"), isolation_level="AUTOCOMMIT")
    try:
        with admin_engine.connect() as conn:
            conn.execute(text(f'CREATE DATABASE "{name}"'))
    except Exception as e:  # pragma: no cover - depends on DB privileges
        pytest.skip(f"cannot create temporary database: {e}")
    url = str(base_url.set(database=name))
    yield url
    create_engine(url).dispose()
    with admin_engine.connect() as conn:
        conn.execute(text(
            f"SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = '{name}' AND pid <> pg_backend_pid()"))
        conn.execute(text(f'DROP DATABASE IF EXISTS "{name}"'))
    admin_engine.dispose()


def _alembic(url):
    cfg = AlembicConfig("alembic.ini")
    cfg.set_main_option("sqlalchemy.url", url)
    return cfg


def _state(url):
    engine = create_engine(url)
    try:
        insp = inspect(engine)
        columns = {c["name"] for c in insp.get_columns("users")} if insp.has_table("users") else set()
        with engine.connect() as conn:
            version = conn.execute(text("SELECT version_num FROM alembic_version")).scalar()
            users = conn.execute(text("SELECT email, must_change_password FROM users")).fetchall() \
                if "must_change_password" in columns else []
        return columns, version, users
    finally:
        engine.dispose()


def _insert_user(url):
    engine = create_engine(url)
    with engine.begin() as conn:
        conn.execute(text(
            "INSERT INTO users (username, email, hashed_password, is_active, is_superuser) "
            "VALUES ('kept', 'kept@example.com', 'x', true, false)"))
    engine.dispose()


def test_a_empty_database(temp_database_url):
    run_migrations(temp_database_url)

    columns, version, _ = _state(temp_database_url)
    assert "must_change_password" in columns
    assert version == HEAD_REVISION


def test_b_create_all_database_without_alembic(temp_database_url):
    # create_all() 시절 스키마 = 최초 revision 스키마 (버전 테이블 없음)
    alembic_command.upgrade(_alembic(temp_database_url), BASELINE_REVISION)
    engine = create_engine(temp_database_url)
    with engine.begin() as conn:
        conn.execute(text("DROP TABLE alembic_version"))
    engine.dispose()
    _insert_user(temp_database_url)

    run_migrations(temp_database_url)

    columns, version, users = _state(temp_database_url)
    assert "must_change_password" in columns
    assert version == HEAD_REVISION
    assert users == [("kept@example.com", False)]


def test_c_alembic_database_at_baseline(temp_database_url):
    alembic_command.upgrade(_alembic(temp_database_url), BASELINE_REVISION)
    _insert_user(temp_database_url)

    run_migrations(temp_database_url)

    columns, version, users = _state(temp_database_url)
    assert version == HEAD_REVISION
    assert users == [("kept@example.com", False)]


def test_d_version_table_without_app_tables(temp_database_url):
    alembic_command.upgrade(_alembic(temp_database_url), BASELINE_REVISION)
    engine = create_engine(temp_database_url)
    with engine.begin() as conn:
        for table in inspect(conn).get_table_names():
            if table != "alembic_version":
                conn.execute(text(f'DROP TABLE "{table}" CASCADE'))
    engine.dispose()

    run_migrations(temp_database_url)

    columns, version, _ = _state(temp_database_url)
    assert "must_change_password" in columns
    assert version == HEAD_REVISION


def test_already_at_head_is_noop(temp_database_url):
    run_migrations(temp_database_url)
    run_migrations(temp_database_url)
    assert _state(temp_database_url)[1] == HEAD_REVISION
