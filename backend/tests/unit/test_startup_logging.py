"""
Unit tests for logger restoration after Alembic migrations (PR-11).

``alembic/env.py`` calls ``fileConfig(alembic.ini)`` with the default
``disable_existing_loggers=True``, which disables every logger that already
exists when ``startup.run_migrations()`` runs. Before the refactor (519649f) the
task-service and datatable-cache loggers were created after migrations and so
stayed enabled; ``restore_loggers_disabled_by_alembic()`` keeps that behavior.
Loggers that were already disabled before the refactor must stay disabled.
"""
import logging
import logging.config
import os

import pytest
from unittest.mock import MagicMock, patch

from app.core import startup
from app.task import service as task_service


ALEMBIC_INI = os.path.join(os.path.dirname(__file__), "..", "..", "alembic.ini")
RESTORED = ("app.task.service", "app.datatable.cache")
STILL_DISABLED = ("app.plugin.service", "app.workflow.service")


@pytest.fixture
def isolated_logging():
    """Snapshot and restore the global logging state touched by fileConfig."""
    root = logging.getLogger()
    root_state = (root.level, list(root.handlers))
    names = RESTORED + STILL_DISABLED + ("sqlalchemy.engine", "alembic")
    saved = {n: (logging.getLogger(n).disabled, logging.getLogger(n).level,
                 list(logging.getLogger(n).handlers)) for n in names}
    yield
    root.setLevel(root_state[0])
    root.handlers[:] = root_state[1]
    for n, (disabled, level, handlers) in saved.items():
        lg = logging.getLogger(n)
        lg.disabled, lg.level = disabled, level
        lg.handlers[:] = handlers


def _apply_alembic_file_config():
    # Same call alembic/env.py makes during run_migrations().
    logging.config.fileConfig(ALEMBIC_INI)


class TestRestoreLoggersDisabledByAlembic:
    def test_file_config_disables_existing_loggers(self, isolated_logging):
        for name in RESTORED + STILL_DISABLED:
            logging.getLogger(name)
        _apply_alembic_file_config()
        for name in RESTORED + STILL_DISABLED:
            assert logging.getLogger(name).disabled is True

    def test_restores_only_loggers_enabled_before_refactor(self, isolated_logging):
        for name in RESTORED + STILL_DISABLED:
            logging.getLogger(name)
        _apply_alembic_file_config()

        startup.restore_loggers_disabled_by_alembic()

        for name in RESTORED:
            assert logging.getLogger(name).disabled is False
        for name in STILL_DISABLED:
            assert logging.getLogger(name).disabled is True

    def test_run_migrations_restores_loggers(self, isolated_logging):
        """run_migrations() re-enables the loggers after every alembic branch."""
        inspector = MagicMock()
        inspector.get_table_names.return_value = ["alembic_version", "users"]

        def fake_upgrade(cfg, revision):
            _apply_alembic_file_config()

        with patch.object(startup, "sa_inspect", return_value=inspector), \
                patch.object(startup.alembic_command, "upgrade", side_effect=fake_upgrade):
            startup.run_migrations()

        for name in RESTORED:
            assert logging.getLogger(name).disabled is False

    def test_task_service_logs_are_emitted(self, isolated_logging):
        _apply_alembic_file_config()
        startup.restore_loggers_disabled_by_alembic()

        # fileConfig replaces the root handlers (incl. pytest's caplog handler),
        # so capture directly on the task-service logger.
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        task_logger = logging.getLogger("app.task.service")
        task_logger.addHandler(handler)
        task_logger.setLevel(logging.INFO)
        try:
            user = MagicMock(id=1, username="tester")
            with patch("app.workflow.compiler.dag_parser.clear_all_caches"):
                task_service.clear_dag_caches(current_user=user)
        finally:
            task_logger.removeHandler(handler)

        assert any("Legacy DAG cache cleared by user 1" in r.getMessage()
                   for r in records)
