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


def _logger_state(lg):
    return lg.disabled, lg.level, lg.propagate, list(lg.handlers)


def _set_logger_state(lg, state):
    lg.disabled, lg.level, lg.propagate, handlers = state
    lg.handlers[:] = handlers


@pytest.fixture
def isolated_logging():
    """Snapshot and restore the whole process-wide logging state.

    fileConfig() disables *every* existing logger (not only the ones these tests
    look at) and replaces the root handlers, so all loggers are snapshotted and
    restored; otherwise later tests (e.g. the ``security`` logger tests) would
    run against disabled loggers. fileConfig() also closes existing handlers;
    file/stream handlers reopen or keep their stream on the next emit, so
    restoring the handler objects is sufficient.
    """
    manager = logging.root.manager
    root = logging.getLogger()
    root_state = _logger_state(root)
    saved = {name: _logger_state(lg) for name, lg in manager.loggerDict.items()
             if isinstance(lg, logging.Logger)}
    yield
    _set_logger_state(root, root_state)
    for name, lg in list(manager.loggerDict.items()):
        if not isinstance(lg, logging.Logger):
            continue
        if name in saved:
            _set_logger_state(lg, saved[name])
        else:
            # Created during the test (e.g. by fileConfig): reset to defaults.
            _set_logger_state(lg, (False, logging.NOTSET, True, []))


def _apply_alembic_file_config(*_args, **_kwargs):
    # Same call alembic/env.py makes during run_migrations().
    logging.config.fileConfig(ALEMBIC_INI)


def _disabled_loggers():
    return {name for name, lg in logging.root.manager.loggerDict.items()
            if isinstance(lg, logging.Logger) and lg.disabled}


class TestRestoreLoggersDisabledByAlembic:
    def test_file_config_disables_existing_loggers(self, isolated_logging):
        for name in RESTORED:
            logging.getLogger(name)
        _apply_alembic_file_config()
        for name in RESTORED:
            assert logging.getLogger(name).disabled is True

    def test_restores_only_loggers_enabled_before_refactor(self, isolated_logging):
        for name in RESTORED:
            logging.getLogger(name)
        _apply_alembic_file_config()
        disabled_by_alembic = _disabled_loggers()
        # The other app loggers (plugin/workflow services, ...) were also disabled
        # before the refactor; they must not be touched.
        assert "app.plugin.service" in disabled_by_alembic

        startup.restore_loggers_disabled_by_alembic()

        assert disabled_by_alembic - _disabled_loggers() == set(RESTORED)

    @pytest.mark.parametrize("tables, alembic_command", [
        (["alembic_version", "users"], "upgrade"),  # alembic-managed DB
        (["users"], "stamp"),                       # legacy create_all() DB
        ([], "upgrade"),                            # fresh install
    ])
    def test_run_migrations_restores_loggers(self, isolated_logging, tables, alembic_command):
        """run_migrations() re-enables the loggers after every alembic branch."""
        inspector = MagicMock()
        inspector.get_table_names.return_value = tables
        for name in RESTORED:
            logging.getLogger(name)

        with patch.object(startup, "sa_inspect", return_value=inspector), \
                patch.object(startup.alembic_command, "upgrade",
                             side_effect=_apply_alembic_file_config) as upgrade, \
                patch.object(startup.alembic_command, "stamp",
                             side_effect=_apply_alembic_file_config) as stamp:
            startup.run_migrations()

        called = {"upgrade": upgrade, "stamp": stamp}
        called[alembic_command].assert_called_once()
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
