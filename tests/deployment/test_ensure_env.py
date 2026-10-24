"""
scripts/ensure-env.sh and compose credential wiring (#146).

Runs the script against a temporary repository root with a stub ``docker`` on PATH,
so no Docker daemon is needed. The compose checks use the real ``docker compose``
(skipped when it is not installed).
"""
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = PROJECT_ROOT / "scripts" / "ensure-env.sh"
CREDENTIAL_KEYS = ["POSTGRES_USER", "POSTGRES_DB", "POSTGRES_PASSWORD",
                   "RABBITMQ_USER", "RABBITMQ_PASSWORD", "REDIS_PASSWORD"]
INSTALL_COMPOSE_FILES = ["docker-compose.cpu.yml", "docker-compose.cpu.amd64.yml", "docker-compose.cpu.arm64.yml",
                         "docker-compose.gpu.yml", "docker-compose.gpu.amd64.yml"]

DOCKER_STUB = """#!/bin/sh
# docker info -> FAKE_DOCKER_UP, docker volume inspect NAME -> NAME in FAKE_VOLUMES
case "$1" in
  info) [ "$FAKE_DOCKER_UP" = "1" ] ;;
  volume)
    for v in $FAKE_VOLUMES; do [ "$v" = "$3" ] && exit 0; done
    exit 1 ;;
  *) exit 1 ;;
esac
"""


def parse_env(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r"^\s*([A-Z_]+)\s*=\s*(.*?)\s*$", line)
        if m:
            values[m.group(1)] = m.group(2)
    return values


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "cellcraft"
    root.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "docker"
    stub.write_text(DOCKER_STUB)
    stub.chmod(0o755)
    return root, bin_dir


def run_script(repo, docker_up=True, volumes="", extra_env=None):
    root, bin_dir = repo
    env = {
        "PATH": f"{bin_dir}:{os.environ['PATH']}",
        "HOME": str(root),
        "CELLCRAFT_ROOT": str(root),
        "FAKE_DOCKER_UP": "1" if docker_up else "0",
        "FAKE_VOLUMES": volumes,
    }
    env.update(extra_env or {})
    return subprocess.run(["bash", str(SCRIPT)], env=env, capture_output=True, text=True, timeout=30)


class TestEnsureEnv:
    def test_new_install_generates_credentials(self, repo):
        root, _ = repo
        result = run_script(repo)
        assert result.returncode == 0, result.stderr

        env_file = root / ".env"
        values = parse_env(env_file)
        assert set(CREDENTIAL_KEYS) <= set(values)
        assert values["POSTGRES_USER"] == "cellcraft_admin" and values["POSTGRES_DB"] == "cellcraft"
        assert values["RABBITMQ_USER"] == "cellcraft"
        secrets = [values[k] for k in ("POSTGRES_PASSWORD", "RABBITMQ_PASSWORD", "REDIS_PASSWORD")]
        assert all(re.fullmatch(r"[0-9a-f]{48}", s) for s in secrets)
        assert len(set(secrets)) == 3
        assert stat.S_IMODE(env_file.stat().st_mode) == 0o600
        # 생성된 값은 출력하지 않는다
        for secret in secrets:
            assert secret not in result.stdout + result.stderr

    def test_rerun_keeps_file_unchanged(self, repo):
        root, _ = repo
        run_script(repo)
        before = (root / ".env").read_bytes()
        result = run_script(repo)
        assert result.returncode == 0
        assert (root / ".env").read_bytes() == before

    def test_existing_values_kept(self, repo):
        root, _ = repo
        existing = ("POSTGRES_USER = custom_user\nPOSTGRES_PASSWORD = \"custom-pass\"\nPOSTGRES_DB = custom_db\n"
                    "RABBITMQ_USER=rmq\nRABBITMQ_PASSWORD='rmq-pass'\nREDIS_PASSWORD=redis-pass\n")
        (root / ".env").write_text(existing)

        result = run_script(repo, docker_up=False)  # 모든 값이 있으면 Docker 없이도 성공
        assert result.returncode == 0, result.stderr
        assert (root / ".env").read_text() == existing

    def test_missing_keys_added_to_partial_file(self, repo):
        # README의 GPU_COUNT 안내처럼 .env가 먼저 만들어진 신규 설치
        root, _ = repo
        (root / ".env").write_text("GPU_COUNT = 4")  # 마지막 줄바꿈 없음
        result = run_script(repo)
        assert result.returncode == 0, result.stderr

        values = parse_env(root / ".env")
        assert values["GPU_COUNT"] == "4"
        assert re.fullmatch(r"[0-9a-f]{48}", values["POSTGRES_PASSWORD"])
        assert set(CREDENTIAL_KEYS) <= set(values)

    def test_legacy_database_volume_uses_previous_account(self, repo):
        result = run_script(repo, volumes="cellcraft_postgres_data")
        assert result.returncode == 0, result.stderr

        values = parse_env(repo[0] / ".env")
        assert values["POSTGRES_USER"] == "cellcraft_admin"
        assert values["POSTGRES_PASSWORD"] == "cellcraft_admin"
        assert values["POSTGRES_DB"] == "cellcraft"
        # RabbitMQ·Redis는 볼륨과 무관하게 새로 생성
        assert values["RABBITMQ_PASSWORD"] != "guest" and values["REDIS_PASSWORD"]
        assert "Existing database volume" in result.stderr

    def test_project_name_normalized_like_compose(self, tmp_path, repo):
        _, bin_dir = repo
        root = tmp_path / "CellCraft.v1"
        root.mkdir()
        result = run_script((root, bin_dir), volumes="cellcraftv1_postgres_data")
        assert result.returncode == 0, result.stderr
        assert parse_env(root / ".env")["POSTGRES_PASSWORD"] == "cellcraft_admin"

    def test_compose_project_name_override(self, repo):
        result = run_script(repo, volumes="other_postgres_data", extra_env={"COMPOSE_PROJECT_NAME": "other"})
        assert result.returncode == 0, result.stderr
        assert parse_env(repo[0] / ".env")["POSTGRES_PASSWORD"] == "cellcraft_admin"

    def test_docker_unavailable_fails_without_writing(self, repo):
        # DB 볼륨 존재 여부를 모르면 비밀번호를 만들지 않는다 (기존 DB와 불일치 방지)
        result = run_script(repo, docker_up=False)
        assert result.returncode == 1
        assert "Docker is not running" in result.stderr
        assert not (repo[0] / ".env").exists()

    def test_env_get_parses_formats(self, repo):
        root, _ = repo
        (root / ".env").write_text('A = spaced  \nB="quoted"\nexport C=exported\r\nA=last\n')
        script = f'source "{SCRIPT}"; for k in A B C D; do printf "%s|" "$(cellcraft_env_get $k)"; done'
        result = subprocess.run(["bash", "-c", script], env={**os.environ, "CELLCRAFT_ROOT": str(root)},
                                capture_output=True, text=True, timeout=10)
        assert result.stdout == "last|quoted|exported||"


@pytest.mark.skipif(shutil.which("docker") is None, reason="docker compose not installed")
class TestComposeCredentials:
    def compose_config(self, compose_file, env_file):
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("POSTGRES_", "RABBITMQ_", "REDIS_", "COMPOSE_"))}
        return subprocess.run(
            ["docker", "compose", "--env-file", str(env_file), "-f", str(PROJECT_ROOT / compose_file),
             "config", "--format", "json"],
            env=env, capture_output=True, text=True, timeout=60,
        )

    @pytest.fixture
    def env_file(self, tmp_path):
        path = tmp_path / ".env"
        path.write_text("POSTGRES_USER=u\nPOSTGRES_PASSWORD=pg-secret\nPOSTGRES_DB=d\n"
                        "RABBITMQ_USER=rmq\nRABBITMQ_PASSWORD=rmq-secret\nREDIS_PASSWORD=redis-secret\n")
        return path

    @pytest.mark.parametrize("compose_file", INSTALL_COMPOSE_FILES)
    def test_internal_services_not_published(self, compose_file, env_file):
        import json
        result = self.compose_config(compose_file, env_file)
        assert result.returncode == 0, result.stderr
        services = json.loads(result.stdout)["services"]

        for name in ("db", "redis", "rabbitmq"):
            assert not services[name].get("ports"), f"{name} publishes ports"
        assert {p["published"] for p in services["backend"]["ports"]} == {"8000"}
        assert {p["published"] for p in services["frontend"]["ports"]} == {"8080"}

        assert services["rabbitmq"]["environment"]["RABBITMQ_DEFAULT_USER"] == "rmq"
        assert services["rabbitmq"]["environment"]["RABBITMQ_DEFAULT_PASS"] == "rmq-secret"
        assert services["redis"]["environment"]["REDIS_PASSWORD"] == "redis-secret"
        for name in ("backend", "celery"):
            env = services[name]["environment"]
            assert env["POSTGRES_PASSWORD"] == "pg-secret"
            assert env["RABBITMQ_USER"] == "rmq" and env["RABBITMQ_PASSWORD"] == "rmq-secret"
            assert env["REDIS_PASSWORD"] == "redis-secret"
        volumes = [v.get("target") for v in services["db"].get("volumes", [])]
        assert "/docker-entrypoint-initdb.d/init.sql" not in volumes

    @pytest.mark.parametrize("compose_file", INSTALL_COMPOSE_FILES + ["docker-compose.dev.yml"])
    def test_missing_env_file_reports_cause(self, compose_file, tmp_path):
        empty = tmp_path / "empty.env"
        empty.write_text("")
        result = self.compose_config(compose_file, empty)
        assert result.returncode != 0
        assert ".env not found" in result.stderr

    def test_dev_binds_internal_ports_to_localhost(self, env_file):
        import json
        result = self.compose_config("docker-compose.dev.yml", env_file)
        assert result.returncode == 0, result.stderr
        services = json.loads(result.stdout)["services"]
        for name in ("redis", "rabbitmq", "test-db"):
            assert all(p.get("host_ip") == "127.0.0.1" for p in services[name]["ports"]), name

    def test_prod_keeps_previous_defaults_without_new_keys(self, tmp_path):
        # 운영 서버의 기존 .env(POSTGRES_*만 있음)로도 이전과 같은 값으로 동작
        import json
        legacy = tmp_path / "legacy.env"
        legacy.write_text("POSTGRES_USER = u\nPOSTGRES_PASSWORD = p\nPOSTGRES_DB = d\n")
        result = self.compose_config("docker-compose.prod.yml", legacy)
        assert result.returncode == 0, result.stderr
        services = json.loads(result.stdout)["services"]
        assert services["rabbitmq"]["environment"]["RABBITMQ_DEFAULT_USER"] == "guest"
        assert services["backend"]["environment"]["RABBITMQ_PASSWORD"] == "guest"
        assert services["redis"]["environment"]["REDIS_PASSWORD"] == ""
