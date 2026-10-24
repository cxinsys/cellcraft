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
from docker.types import Mount

from app.core.config import settings
from app.workflow.compiler import snakefile as sf


HOST_BACKEND = "/srv/cellcraft/backend"


def _bind(target, source, read_only):
    return Mount(target=target, source=source, type="bind", read_only=read_only)


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


class TestToContainerPath:
    @pytest.mark.parametrize("path", [
        "./user/john/workflow_1/algorithm_2/Snakefile",
        "user/john/workflow_1/algorithm_2/Snakefile",
    ])
    def test_with_or_without_dot_prefix(self, path):
        assert sf.to_container_path(path) == "/workspace/user/john/workflow_1/algorithm_2/Snakefile"


class TestBuildPluginMounts:
    def test_mounts_only_user_dir_and_shared_inputs(self, tmp_path):
        (tmp_path / "tutorials").mkdir()
        mounts = sf.build_plugin_mounts(HOST_BACKEND, tmp_path, "user/john")
        assert mounts == [
            _bind("/workspace/user/john", f"{HOST_BACKEND}/user/john", read_only=False),
            _bind("/workspace/tutorials", f"{HOST_BACKEND}/tutorials", read_only=True),
        ]

    def test_missing_shared_dir_not_mounted(self, tmp_path):
        mounts = sf.build_plugin_mounts(HOST_BACKEND, tmp_path, "user/john")
        assert [m["Target"] for m in mounts] == ["/workspace/user/john"]

    def test_accepts_path_objects(self, tmp_path):
        mounts = sf.build_plugin_mounts(tmp_path, tmp_path, "user/john")
        assert mounts[0]["Source"] == f"{tmp_path}/user/john"

    def test_colon_in_username_kept_intact(self, tmp_path):
        # volumes 문자열(src:dst:mode)과 달리 Mount는 ':'를 구분자로 해석하지 않는다
        mounts = sf.build_plugin_mounts(HOST_BACKEND, tmp_path, "user/a:b")
        assert mounts[0]["Source"] == f"{HOST_BACKEND}/user/a:b"
        assert mounts[0]["Target"] == "/workspace/user/a:b"


class TestSettingsDefaults:
    def test_network_default_is_bridge(self):
        assert type(settings).__fields__["PLUGIN_CONTAINER_NETWORK"].default == "bridge"

    def test_hardening_default_enabled(self):
        assert type(settings).__fields__["PLUGIN_CONTAINER_HARDENING"].default is True


class TestSecurityOptions:
    def test_hardening_enabled(self, monkeypatch):
        monkeypatch.setattr(settings, "PLUGIN_CONTAINER_HARDENING", True)
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


def _run_exec_in_plugin(tmp_path, monkeypatch, client, snakefile_path, task_id=None):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOSTNAME", "celery-container")
    monkeypatch.setattr(settings, "PLUGIN_CONTAINER_NETWORK", "bridge")
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
        return sf.exec_in_plugin("MyPlugin", snakefile_path, ["out.txt"], str(tmp_path), task_id)


class TestExecInPluginContainerConfig:
    SNAKEFILE = "./user/john/workflow_1/algorithm_2/Snakefile"

    @pytest.fixture(autouse=True)
    def _workspace(self, tmp_path):
        (tmp_path / "user/john/workflow_1/algorithm_2").mkdir(parents=True)

    def test_container_isolated_from_backend_and_celery(self, tmp_path, monkeypatch, docker_client):
        monkeypatch.setattr(settings, "PLUGIN_CONTAINER_HARDENING", True)
        result = _run_exec_in_plugin(tmp_path, monkeypatch, docker_client, self.SNAKEFILE)

        assert result["returncode"] == 0
        config = docker_client.containers.run.call_args.kwargs
        assert "volumes_from" not in config
        assert "volumes" not in config
        assert config["mounts"] == [
            _bind("/workspace/user/john", f"{HOST_BACKEND}/user/john", read_only=False),
            _bind("/workspace/tutorials", f"{HOST_BACKEND}/tutorials", read_only=True),
        ]
        assert config["working_dir"] == "/workspace"
        assert config["network"] == "bridge"
        assert config["cap_drop"] == ["ALL"]
        assert config["security_opt"] == ["no-new-privileges"]

    def test_hardening_disabled_keeps_isolation(self, tmp_path, monkeypatch, docker_client):
        monkeypatch.setattr(settings, "PLUGIN_CONTAINER_HARDENING", False)
        _run_exec_in_plugin(tmp_path, monkeypatch, docker_client, self.SNAKEFILE)

        config = docker_client.containers.run.call_args.kwargs
        assert "cap_drop" not in config
        assert "security_opt" not in config
        assert "volumes_from" not in config
        assert [m["Target"] for m in config["mounts"]] == ["/workspace/user/john", "/workspace/tutorials"]

    def test_snakefile_path_unchanged_inside_container(self, tmp_path, monkeypatch, docker_client):
        _run_exec_in_plugin(tmp_path, monkeypatch, docker_client, self.SNAKEFILE)

        bash_cmd = docker_client.api.exec_create.call_args.args[1][2]
        assert "--snakefile /workspace/user/john/workflow_1/algorithm_2/Snakefile" in bash_cmd

    def test_task_logs_stay_inside_user_dir(self, tmp_path, monkeypatch, docker_client):
        result = _run_exec_in_plugin(
            tmp_path, monkeypatch, docker_client, self.SNAKEFILE, task_id="task-1"
        )

        assert result["returncode"] == 0
        node_dir = tmp_path / "user/john/workflow_1/algorithm_2"
        assert (node_dir / "logs").is_symlink()
        assert (node_dir / "logs").resolve() == (node_dir / "executions/task-1/logs").resolve()
        assert (node_dir / "executions/task-1/meta.yml").exists()


class TestExecInPluginRejectsOutsideUserDir:
    def test_rejected_before_any_side_effect(self, tmp_path, monkeypatch, docker_client):
        result = _run_exec_in_plugin(
            tmp_path, monkeypatch, docker_client, "./user/../workflow_1/algorithm_2/Snakefile",
            task_id="task-1",
        )

        assert result["returncode"] == 1
        docker_client.containers.run.assert_not_called()
        # 로그 디렉터리·symlink 등 파일 시스템 변경이 일어나지 않아야 한다
        assert not (tmp_path / "workflow_1").exists()
