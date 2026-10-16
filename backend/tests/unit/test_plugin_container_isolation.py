"""
Unit tests for plugin execution container isolation (#136).

Pins how ``exec_in_plugin`` builds the plugin container: only the executing
user's directory (rw) and shared inputs (ro) are mounted, the Celery container's
volumes are not inherited, and the container runs on the configured network
with dropped capabilities. Docker itself is fully mocked.
"""
import contextlib
from unittest.mock import MagicMock, patch

import pytest

from app.core.config import settings
from app.workflow.compiler import snakefile as sf


HOST_BACKEND = "/srv/cellcraft/backend"


class TestResolveUserDir:
    @pytest.mark.parametrize("path, expected", [
        ("./user/john/workflow_1/algorithm_2/Snakefile", "user/john"),
        ("user/john/workflow_1/visualization_3/Snakefile", "user/john"),
        ("./user/john/./workflow_1/algorithm_2/Snakefile", "user/john"),
    ])
    def test_valid_paths(self, path, expected):
        assert sf.resolve_user_dir(path) == expected

    @pytest.mark.parametrize("path", [
        "./user/../plugin/official/TENET/Snakefile",
        "./user/john/../../app/Snakefile",
        "/app/user/john/workflow_1/Snakefile",
        "./tutorials/Snakefile",
        "./user/Snakefile",
        "Snakefile",
    ])
    def test_paths_outside_user_dir_rejected(self, path):
        with pytest.raises(ValueError):
            sf.resolve_user_dir(path)


class TestBuildPluginVolumes:
    def test_mounts_only_user_dir_and_shared_inputs(self, tmp_path):
        (tmp_path / "tutorials").mkdir()
        volumes = sf.build_plugin_volumes(
            HOST_BACKEND, tmp_path, "./user/john/workflow_1/algorithm_2/Snakefile"
        )
        assert volumes == {
            f"{HOST_BACKEND}/user/john": {"bind": "/workspace/user/john", "mode": "rw"},
            f"{HOST_BACKEND}/tutorials": {"bind": "/workspace/tutorials", "mode": "ro"},
        }
        assert HOST_BACKEND not in volumes

    def test_missing_shared_dir_not_mounted(self, tmp_path):
        volumes = sf.build_plugin_volumes(
            HOST_BACKEND, tmp_path, "./user/john/workflow_1/algorithm_2/Snakefile"
        )
        assert list(volumes) == [f"{HOST_BACKEND}/user/john"]

    def test_accepts_path_objects(self, tmp_path):
        volumes = sf.build_plugin_volumes(
            tmp_path, tmp_path, "./user/john/workflow_1/algorithm_2/Snakefile"
        )
        assert f"{tmp_path}/user/john" in volumes


class TestSecurityOptions:
    def test_hardening_enabled_by_default(self):
        assert sf.plugin_container_security_options() == {
            "cap_drop": ["ALL"],
            "security_opt": ["no-new-privileges"],
        }

    def test_hardening_can_be_disabled(self, monkeypatch):
        monkeypatch.setattr(settings, "PLUGIN_CONTAINER_HARDENING", False)
        assert sf.plugin_container_security_options() == {}


@pytest.fixture
def docker_client():
    """Mocked docker client whose celery container mounts the host backend at /app."""
    client = MagicMock()
    celery_container = MagicMock()
    celery_container.attrs = {"Mounts": [
        {"Destination": "/var/run/docker.sock", "Source": "/var/run/docker.sock"},
        {"Destination": "/app", "Source": HOST_BACKEND},
    ]}
    client.containers.get.return_value = celery_container
    plugin_container = MagicMock()
    plugin_container.id = "plugin-container-id"
    plugin_container.status = "running"
    client.containers.run.return_value = plugin_container
    client.api.exec_create.return_value = {"Id": "exec-id"}
    client.api.exec_start.return_value = iter([b"done"])
    client.api.exec_inspect.return_value = {"ExitCode": 0}
    return client


def _run_exec_in_plugin(tmp_path, monkeypatch, client, snakefile_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOSTNAME", "celery-container")
    (tmp_path / "tutorials").mkdir(exist_ok=True)

    db = MagicMock()
    db.query.return_value.filter_by.return_value.first.return_value = None  # local plugin

    @contextlib.contextmanager
    def fake_session():
        yield db

    with patch.object(sf.docker, "from_env", return_value=client), \
         patch("app.db.session.get_db_session", fake_session), \
         patch.object(sf, "wait_for_container_ready"), \
         patch.object(sf.container_manager, "register_container"), \
         patch.object(sf.container_manager, "unregister_container"), \
         patch.object(sf.container_manager, "_is_cleanup_in_progress", return_value=False):
        return sf.exec_in_plugin("MyPlugin", snakefile_path, ["out.txt"], str(tmp_path))


class TestExecInPluginContainerConfig:
    def test_container_isolated_from_backend_and_celery(self, tmp_path, monkeypatch, docker_client):
        snakefile = "./user/john/workflow_1/algorithm_2/Snakefile"
        (tmp_path / "user/john/workflow_1/algorithm_2").mkdir(parents=True)

        result = _run_exec_in_plugin(tmp_path, monkeypatch, docker_client, snakefile)

        assert result["returncode"] == 0
        config = docker_client.containers.run.call_args.kwargs
        assert "volumes_from" not in config
        assert config["volumes"] == {
            f"{HOST_BACKEND}/user/john": {"bind": "/workspace/user/john", "mode": "rw"},
            f"{HOST_BACKEND}/tutorials": {"bind": "/workspace/tutorials", "mode": "ro"},
        }
        assert config["working_dir"] == "/workspace"
        assert config["network"] == settings.PLUGIN_CONTAINER_NETWORK
        assert config["cap_drop"] == ["ALL"]
        assert config["security_opt"] == ["no-new-privileges"]

    def test_snakefile_path_unchanged_inside_container(self, tmp_path, monkeypatch, docker_client):
        snakefile = "./user/john/workflow_1/algorithm_2/Snakefile"
        (tmp_path / "user/john/workflow_1/algorithm_2").mkdir(parents=True)

        _run_exec_in_plugin(tmp_path, monkeypatch, docker_client, snakefile)

        bash_cmd = docker_client.api.exec_create.call_args.args[1][2]
        assert "--snakefile /workspace/user/john/workflow_1/algorithm_2/Snakefile" in bash_cmd

    def test_snakefile_outside_user_dir_not_executed(self, tmp_path, monkeypatch, docker_client):
        result = _run_exec_in_plugin(
            tmp_path, monkeypatch, docker_client, "./user/../plugin/Snakefile"
        )

        assert result["returncode"] == 1
        docker_client.containers.run.assert_not_called()
