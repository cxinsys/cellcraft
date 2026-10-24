"""Application startup logic for the FastAPI web layer.

Holds the pieces that ran at import time / on the ``startup`` event in the old
``main.py``: DB migrations, signal handlers, the plugin-system startup routine,
and Docker image pulling. Logic is moved verbatim from ``main.py`` — behavior
(including print-based logging) is unchanged.
"""
import signal
import sys
import os
import asyncio
from concurrent.futures import ThreadPoolExecutor
import logging

import docker
from alembic.config import Config as AlembicConfig
from alembic import command as alembic_command
from sqlalchemy import inspect as sa_inspect

from app.core.config import settings
from app.shared.docker import container_manager
from app.plugin.sync import PluginSyncManager, initialize_plugins_from_csv
from app.plugin.versioning import PluginVersionValidator
from app import models
from app.db.session import engine, get_db_session

# Docker pull 작업을 위한 ThreadPoolExecutor (메인 이벤트 루프 블로킹 방지)
_docker_executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix="docker-pull")

# alembic/env.py의 fileConfig(alembic.ini)는 기본값 disable_existing_loggers=True로
# 그 시점에 이미 생성된 logger를 모두 비활성화한다. 리팩토링 전(519649f)에는 아래
# logger들이 마이그레이션 이후(요청 시점·지연 import 시점)에 생성되어 출력되었으나,
# 현재는 모듈 최상단 import로 마이그레이션 전에 생성된다. 로그 출력 범위를 리팩토링 전과
# 동일하게 유지하기 위해 이 logger들만 다시 활성화한다.
_LOGGERS_ENABLED_BEFORE_REFACTOR = ("app.task.service", "app.datatable.cache")


def restore_loggers_disabled_by_alembic():
    for name in _LOGGERS_ENABLED_BEFORE_REFACTOR:
        logging.getLogger(name).disabled = False


# create_all() 시절 스키마와 동일한 최초 revision (상태 B DB의 기준점)
BASELINE_REVISION = "0001_initial"


def run_migrations(database_url: str = None):
    """Alembic 마이그레이션을 실행하여 DB 스키마를 최신 상태로 유지한다.

    기존 create_all() 대신 Alembic을 Single Source of Truth로 사용.
    - 첫 설치 (상태 A): 빈 DB → upgrade head → 모든 테이블 생성
    - 기존 DB (상태 B): create_all()로 생성된 DB(= 최초 revision 스키마) → stamp baseline → upgrade head
    - Alembic DB (상태 C): 이미 관리 중 → upgrade head → 미적용 revision만 적용

    Args:
        database_url: 대상 DB (기본값: 앱이 사용하는 engine의 DB, 테스트에서 별도 DB 지정용)
    """
    from app.db.base_registry import Base

    # 상태 판별(inspector)과 마이그레이션(alembic)은 항상 같은 DB를 대상으로 한다
    target_engine = engine
    if database_url and database_url != engine.url.render_as_string(hide_password=False):
        from sqlalchemy import create_engine
        target_engine = create_engine(database_url)
    database_url = target_engine.url.render_as_string(hide_password=False)
    alembic_cfg = AlembicConfig("alembic.ini")
    alembic_cfg.set_main_option("sqlalchemy.url", database_url.replace("%", "%%"))

    try:
        existing_tables = set(sa_inspect(target_engine).get_table_names())
    finally:
        if target_engine is not engine:
            target_engine.dispose()

    has_alembic_version = "alembic_version" in existing_tables
    has_app_tables = "users" in existing_tables
    present_app_tables = existing_tables & set(Base.metadata.tables)

    if has_alembic_version and not present_app_tables:
        # 버전 기록만 남고 앱 테이블이 하나도 없는 DB (테스트 정리·수동 삭제) → 빈 DB로 처리
        print("DB: 앱 테이블 없는 Alembic DB 감지 → 버전 초기화 후 upgrade head", flush=True)
        alembic_command.stamp(alembic_cfg, "base")
        alembic_command.upgrade(alembic_cfg, "head")
    elif not has_app_tables and present_app_tables:
        # 일부 앱 테이블만 남은 비정상 DB → 자동 처리하지 않고 중단 (버전 기록 보존)
        raise RuntimeError(
            "Database schema is incomplete: 'users' table is missing while other application tables "
            f"exist ({', '.join(sorted(present_app_tables))}). Restore the database or recreate it."
        )
    elif has_alembic_version:
        # 상태 C: 이미 Alembic으로 관리 중인 DB → 미적용 revision만 적용
        print("DB: Alembic 관리 DB 감지 → upgrade head", flush=True)
        alembic_command.upgrade(alembic_cfg, "head")
    elif has_app_tables:
        # 상태 B: create_all()로 생성된 기존 DB → 최초 revision으로 표시 후 이후 revision 적용
        print(f"DB: 기존 create_all() DB 감지 → stamp {BASELINE_REVISION} → upgrade head", flush=True)
        alembic_command.stamp(alembic_cfg, BASELINE_REVISION)
        alembic_command.upgrade(alembic_cfg, "head")
    else:
        # 상태 A: 완전히 새로운 설치 → 전체 마이그레이션 실행
        print("DB: 빈 DB 감지 → upgrade head (전체 마이그레이션)", flush=True)
        alembic_command.upgrade(alembic_cfg, "head")

    restore_loggers_disabled_by_alembic()


def setup_signal_handlers():
    """전역 시그널 핸들러 설정"""
    def signal_handler(signum, frame):
        print(f"Received signal {signum}. Cleaning up all containers...")
        try:
            container_manager.cleanup_all_task_containers()
        except Exception as e:
            print(f"Error during container cleanup: {e}")
        finally:
            sys.exit(0)

    # 시그널 핸들러 등록
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)


# 서버 시작 시 플러그인 초기화 및 이미지 Pull
async def on_startup():
    logger = logging.getLogger(__name__)

    print("=" * 80, flush=True)
    print("CELLCRAFT SERVER STARTUP", flush=True)
    print("=" * 80, flush=True)

    # 초기 관리자 계정 (superuser가 없을 때만 생성, 기본 비밀번호 사용 계정 표시)
    print("\n0. Admin Account Initialization...", flush=True)
    try:
        from app.core.bootstrap import bootstrap_admin
        if os.environ.get("TESTING") == "1":
            print("Skipped in test environment", flush=True)
        else:
            with get_db_session() as db:
                bootstrap_admin(db)
    except Exception as e:
        logger.error(f"Admin account initialization failed: {e}", exc_info=True)
        print(f"Admin account initialization failed: {e}", flush=True)

    # Plugin synchronization and consistency check
    print("\n1. Plugin System Initialization...", flush=True)
    try:
        sync_manager = PluginSyncManager()
        validator = PluginVersionValidator()

        # Get version information from file (no git required)
        print("   Reading plugin version information...", flush=True)
        version_status = sync_manager.get_sync_status()

        if "error" in version_status:
            print(f"   ⚠️  Warning: Could not read version info: {version_status['error']}", flush=True)
            print("   Using default version settings", flush=True)
        else:
            version_info = version_status.get("version_info", {})
            print(f"   Plugin version: {version_status.get('repository_version', 'unknown')}", flush=True)
            if version_info:
                print(f"   Build time: {version_info.get('build_time', 'unknown')}", flush=True)
                print(f"   Commit: {version_info.get('commit', 'unknown')}", flush=True)

        # Perform consistency check (simplified - skip Docker registry checks)
        print("\n2. Version Consistency Check...", flush=True)
        try:
            consistency_result = validator.validate_consistency()

            if consistency_result.get("consistent", False):
                print("   ✅ All components are in sync", flush=True)
            else:
                print("   ⚠️  Version differences detected (will be handled during sync)", flush=True)
                # Don't block on version issues - let sync handle it

                # Attempt automatic synchronization
                print("\n3. Automatic Synchronization...", flush=True)
                try:
                    sync_result = sync_manager.sync_plugins_to_database()

                    if sync_result.get("success", False):
                        print("   ✅ Database synchronized successfully", flush=True)
                        print(f"   Updated to version: {sync_result.get('version', 'unknown')}", flush=True)
                    else:
                        print("   ⚠️  Synchronization incomplete (non-critical)", flush=True)
                        print("   Using existing database entries", flush=True)
                except Exception as sync_e:
                    print(f"   ⚠️  Sync skipped: {sync_e}", flush=True)
                    print("   Using existing database entries", flush=True)

        except Exception as e:
            print(f"   ⚠️  Consistency check skipped (non-critical): {e}", flush=True)
            logger.warning(f"Version consistency check had issues: {e}")

        # Ensure plugins are initialized
        if not version_status.get("database_plugin_count", 0):
            print("   Initializing plugins from CSV...", flush=True)
            initialized_count = initialize_plugins_from_csv("./plugin/official/plugins.csv")
            print(f"   ✅ {initialized_count} plugins initialized", flush=True)
        else:
            print(f"   ✅ {version_status.get('database_plugin_count', 0)} plugins already in database", flush=True)

    except Exception as e:
        print(f"   ⚠️  Plugin initialization warning: {e}", flush=True)
        logger.warning(f"Plugin system initialization had issues: {e}")
        # Always try to initialize from CSV as fallback
        try:
            print("   Attempting fallback initialization from plugins.csv...", flush=True)
            fallback_count = initialize_plugins_from_csv("./plugin/official/plugins.csv")
            print(f"   ✅ Fallback initialization successful: {fallback_count} plugins", flush=True)
        except Exception as fallback_error:
            print(f"   ❌ Critical: Could not initialize plugins: {fallback_error}", flush=True)
            logger.error(f"Failed to initialize plugins: {fallback_error}")

    # Docker image pulling - always attempt this
    print("\n4. Docker Images Check...", flush=True)
    try:
        await check_and_pull_official_plugin_images()
        print("   Docker image check completed", flush=True)
    except Exception as e:
        print(f"   ⚠️  Docker image check had issues (non-critical): {e}", flush=True)
        logger.warning(f"Docker image check had issues: {e}")

    print("\n" + "=" * 80, flush=True)
    print("STARTUP COMPLETE", flush=True)
    print("=" * 80, flush=True)


def _sync_pull_image(client, image_uri: str) -> tuple[bool, str]:
    """동기적으로 Docker 이미지를 pull (ThreadPoolExecutor에서 실행, 진행률 표시)"""
    try:
        # 스트리밍으로 진행률 표시
        layer_progress = {}  # layer_id -> (current, total)
        last_status = ""

        for chunk in client.api.pull(image_uri, stream=True, decode=True):
            status = chunk.get("status", "")
            layer_id = chunk.get("id", "")
            progress_detail = chunk.get("progressDetail", {})

            # 다운로드/추출 진행률 추적
            if layer_id and progress_detail:
                current = progress_detail.get("current", 0)
                total = progress_detail.get("total", 0)
                if total > 0:
                    layer_progress[layer_id] = (current, total)

            # 전체 진행률 계산 및 출력 (Downloading/Extracting 상태일 때만)
            if status in ("Downloading", "Extracting") and layer_progress:
                total_current = sum(p[0] for p in layer_progress.values())
                total_size = sum(p[1] for p in layer_progress.values())
                if total_size > 0:
                    percent = (total_current / total_size) * 100
                    mb_current = total_current / (1024 * 1024)
                    mb_total = total_size / (1024 * 1024)
                    progress_line = f"      {status}: {percent:5.1f}% ({mb_current:.1f}/{mb_total:.1f} MB)"
                    # 같은 줄에 덮어쓰기 (carriage return)
                    print(f"\r{progress_line}", end="", flush=True)
                    last_status = status

            # 단계 변경 시 줄바꿈
            elif status and status != last_status and status not in ("Downloading", "Extracting"):
                if last_status in ("Downloading", "Extracting"):
                    print(flush=True)  # 진행률 줄 마무리
                last_status = status

        # 완료 후 줄바꿈
        if last_status in ("Downloading", "Extracting"):
            print(flush=True)

        return (True, f"   ✓ Successfully pulled {image_uri}")
    except docker.errors.APIError as e:
        error_msg = str(e).lower()
        if "not found" in error_msg or "404" in error_msg:
            return (False, f"   ⚠ Image {image_uri} not found in registry")
        elif "unauthorized" in error_msg or "401" in error_msg:
            return (False, f"   ⚠ Cannot access {image_uri} (authentication may be required)")
        else:
            return (False, f"   ⚠ Failed to pull {image_uri}: {e}")
    except Exception as e:
        return (False, f"   ⚠ Could not pull {image_uri}: {e}")


async def check_and_pull_official_plugin_images():
    """Check and pull official plugin Docker images (non-blocking)"""
    client = None  # Docker 클라이언트 누수 방지를 위해 초기화
    loop = asyncio.get_event_loop()

    try:
        with get_db_session() as db:
            client = docker.from_env()

            # Get all official plugins from database
            plugins = db.query(models.Plugin).filter_by(source="official").all()

            # Filter out GPU-only plugins when in CPU-only mode
            cpu_only = os.getenv("CPU_ONLY", "false").lower() == "true"
            if cpu_only:
                gpu_only_plugins = {"FastSCODE", "FastTENET"}
                filtered_plugins = [p for p in plugins if p.name not in gpu_only_plugins]
                if len(filtered_plugins) != len(plugins):
                    print(f"   CPU-only mode: Filtered out {len(plugins) - len(filtered_plugins)} GPU-only plugins", flush=True)
                plugins = filtered_plugins
                print(f"   Found {len(plugins)} CPU-compatible plugins to check", flush=True)
            else:
                print(f"   Found {len(plugins)} official plugins to check", flush=True)

            for plugin in plugins:
                try:
                    # Generate simple image URI without registry client
                    plugin_name_lower = plugin.name.lower()
                    version = plugin.version or "latest"
                    image_uri = f"ghcr.io/cxinsys/cellcraft-{plugin_name_lower}:{version}"

                    # Check if image exists locally (빠른 동기 호출)
                    try:
                        client.images.get(image_uri)
                        print(f"   ✓ Image {image_uri} already exists locally", flush=True)
                    except docker.errors.ImageNotFound:
                        # Try to pull from registry (비동기로 실행하여 이벤트 루프 블로킹 방지)
                        print(f"   ⬇ Pulling {image_uri}...", flush=True)

                        # ThreadPoolExecutor에서 동기 pull 실행
                        _, message = await loop.run_in_executor(
                            _docker_executor,
                            _sync_pull_image,
                            client,
                            image_uri
                        )
                        print(message, flush=True)

                except Exception as e:
                    print(f"   ⚠ Error processing plugin {plugin.name}: {e}", flush=True)
                    continue

    except docker.errors.DockerException as e:
        print(f"   ⚠ Docker not available: {e}", flush=True)
    except Exception as e:
        print(f"   ⚠ Error checking/pulling images: {e}", flush=True)
    finally:
        # Docker 클라이언트 연결 해제 - TCP 소켓 누수 방지
        if client:
            try:
                client.close()
            except Exception:
                pass
