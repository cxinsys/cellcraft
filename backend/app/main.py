from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware
import signal
import sys
import os

from app.routes.api import api_router
from app.common.config import settings
from app.common.utils.celery_utils import create_celery
from app.common.utils.docker_utils import container_manager
from app.database import models
from app.database.conn import engine, initialize_plugins_from_csv, get_new_engine_and_session
from celery.signals import worker_shutting_down
import docker
from app.common.utils.github_registry_client import GitHubRegistryClient
from app.common.utils.plugin_sync_manager import PluginSyncManager
from app.common.utils.plugin_version_validator import PluginVersionValidator
import logging

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

# Celery 이벤트 핸들러를 정의합니다
def on_worker_shut_down(sender=None, conf=None, **kwargs):
    print("Worker is shutting down. Cleaning up containers and trying to reconnect...")
    try:
        # 컨테이너 정리
        container_manager.cleanup_all_task_containers()
        
        # 재연결 시도
        with celery.connection() as connection:
            connection.ensure_connection(max_retries=3)
            print("재연결에 성공했습니다.")
    except Exception as e:
        print(f"Worker shutdown 처리 중 오류 발생: {e}")

models.Base.metadata.create_all(bind=engine)
global_engine = engine

app = FastAPI(
    title=settings.PROJECT_NAME
)

app.add_middleware(
    CORSMiddleware,
    # allow_origins=[str(origin) for origin in settings.BACKEND_CORS_ORIGINS],
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.celery_app = create_celery()
# Celery 이벤트 핸들러를 Celery 애플리케이션 인스턴스에 연결합니다
worker_shutting_down.connect(on_worker_shut_down, sender=app.celery_app)

app.include_router(api_router, prefix=settings.ROUTES_STR)

celery = app.celery_app

# 시그널 핸들러 설정
setup_signal_handlers()

def get_celery_app():
    return celery

# 서버 시작 시 플러그인 초기화 및 이미지 Pull
@app.on_event("startup")
async def startup_event():
    logger = logging.getLogger(__name__)
    
    print("=" * 80)
    print("CELLCRAFT SERVER STARTUP")
    print("=" * 80)
    
    # Plugin synchronization and consistency check
    print("\n1. Plugin System Initialization...")
    try:
        sync_manager = PluginSyncManager()
        validator = PluginVersionValidator()
        
        # Get version information from file (no git required)
        print("   Reading plugin version information...")
        version_status = sync_manager.get_sync_status()
        
        if "error" in version_status:
            print(f"   ⚠️  Warning: Could not read version info: {version_status['error']}")
            print("   Using default version settings")
        else:
            version_info = version_status.get("version_info", {})
            print(f"   Plugin version: {version_status.get('repository_version', 'unknown')}")
            if version_info:
                print(f"   Build time: {version_info.get('build_time', 'unknown')}")
                print(f"   Commit: {version_info.get('commit', 'unknown')}")
        
        # Perform consistency check
        print("\n2. Version Consistency Check...")
        try:
            consistency_result = validator.validate_consistency()
            
            if consistency_result.get("consistent", False):
                print("   ✅ All components are in sync")
            else:
                print("   ⚠️  Version inconsistencies detected")
                issues = consistency_result.get("issues", [])
                for issue in issues[:3]:  # Show first 3 issues
                    print(f"      - {issue}")
                
                if len(issues) > 3:
                    print(f"      ... and {len(issues) - 3} more issues")
                
                # Attempt automatic synchronization
                print("\n3. Automatic Synchronization...")
                sync_result = sync_manager.sync_plugins_to_database()
                
                if sync_result.get("success", False):
                    print("   ✅ Database synchronized successfully")
                    print(f"   Updated to version: {sync_result.get('version', 'unknown')}")
                else:
                    print("   ⚠️  Synchronization incomplete: {sync_result.get('error', 'Unknown error')}")
                    print("   Using existing database entries")
        
        except Exception as e:
            print(f"   ⚠️  Consistency check skipped: {e}")
            logger.warning(f"Version consistency check failed: {e}")
        
        # Ensure plugins are initialized
        if not version_status.get("database_plugin_count", 0):
            print("   Initializing plugins from CSV...")
            initialize_plugins_from_csv("./plugin/official/plugins.csv")
            print("   ✅ Plugins initialized")
            
    except Exception as e:
        print(f"   ⚠️  Plugin initialization warning: {e}")
        logger.warning(f"Plugin system initialization had issues: {e}")
        # Always try to initialize from CSV as fallback
        try:
            print("   Attempting fallback initialization from plugins.csv...")
            initialize_plugins_from_csv("./plugin/official/plugins.csv")
            print("   ✅ Fallback initialization successful")
        except Exception as fallback_error:
            print(f"   ❌ Critical: Could not initialize plugins: {fallback_error}")
            logger.error(f"Failed to initialize plugins: {fallback_error}")
    
    # Docker image pulling
    print("\n4. Docker Images Check...")
    await check_and_pull_official_plugin_images()
    
    print("\n" + "=" * 80)
    print("STARTUP COMPLETE")
    print("=" * 80)

async def check_and_pull_official_plugin_images():
    """Check and pull official plugin Docker images from GitHub Registry"""
    db = get_new_engine_and_session()
    
    try:
        client = docker.from_env()
        registry = GitHubRegistryClient()
        
        # Get all official plugins from database
        plugins = db.query(models.Plugin).filter_by(source="official").all()
        print(f"Found {len(plugins)} official plugins to check")
        
        for plugin in plugins:
            try:
                # Generate image URI
                image_uri = registry.get_image_uri(
                    plugin.name.lower(), 
                    plugin.version or "latest"
                )
                
                # Check if image exists locally
                try:
                    client.images.get(image_uri)
                    print(f"✓ Image {image_uri} already exists locally")
                except docker.errors.ImageNotFound:
                    # Pull from registry
                    print(f"⬇ Pulling {image_uri}...")
                    try:
                        client.images.pull(image_uri)
                        print(f"✓ Successfully pulled {image_uri}")
                    except docker.errors.APIError as e:
                        if "not found" in str(e).lower():
                            print(f"⚠ Image {image_uri} not found in registry, will fallback to local build if needed")
                        else:
                            print(f"❌ Failed to pull {image_uri}: {e}")
                    except Exception as e:
                        print(f"❌ Unexpected error pulling {image_uri}: {e}")
                
            except Exception as e:
                print(f"❌ Error processing plugin {plugin.name}: {e}")
                continue
                
    except docker.errors.DockerException as e:
        print(f"⚠ Docker not available: {e}")
    except Exception as e:
        print(f"❌ Error checking/pulling images: {e}")
    finally:
        db.close()