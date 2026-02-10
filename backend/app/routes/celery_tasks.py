from celery import shared_task, Task
from datetime import datetime, timezone
from celery.worker.request import Request
import os
from pathlib import Path
import time
from typing import List
from billiard import Pool, cpu_count
import logging

from app.common.utils.snakemake_utils import snakemakeProcess
from app.common.utils.docker_utils import container_manager
from app.common.utils import plugin_utils
from app.database.crud.crud_task import start_task, end_task, record_plugin_image_uri
from app.common.utils.cache_utils import save_result_to_cache

logger = logging.getLogger('celery.custom')

class MyRequest(Request):
    """Custom request to detect RuntimeError and ensure task failure is logged correctly."""

    def on_failure(self, exc_info, send_failed_event=True, return_ok=False):
        """Handle task failure and log necessary information."""
        super().on_failure(
            exc_info,
            send_failed_event=send_failed_event,
            return_ok=return_ok
        )

        exception = exc_info.exception
        if isinstance(exception, RuntimeError):
            logger.warning(f"RuntimeError detected in task {self.task.name}: {exception}")
        else:
            logger.warning(f"Failure detected in task {self.task.name}: {exception}")

class MyTask(Task):
    Request = MyRequest  # Custom Request class 적용

    def before_start(self, task_id, args, kwargs):
        start_time = datetime.now(timezone.utc)
        print(f'Task {task_id} started at {start_time}')
        user_id = kwargs.get('user_id')
        workflow_id = kwargs.get('workflow_id')
        algorithm_id = kwargs.get('algorithm_id')
        plugin_name = kwargs.get('plugin_name')
        task_type = kwargs.get('task_type')
        
        # Generate plugin_image_uri and get plugin_id if plugin_name is provided
        plugin_image_uri = None
        plugin_id = None
        if plugin_name:
            try:
                from app.common.utils.plugin_utils import generate_plugin_image_uri
                from app.database.conn import get_db_session
                from app.database import models

                # Get plugin information from database
                with get_db_session() as db:
                    plugin = db.query(models.Plugin).filter_by(name=plugin_name).first()
                    if plugin:
                        source = plugin.source
                        version = plugin.version
                        plugin_id = plugin.id
                        print(f'Found plugin in database: {plugin_name} (id: {plugin_id}, source: {source})')
                    else:
                        source = "local"  # Default to local if not found
                        version = None
                        print(f'Plugin {plugin_name} not found in database, using defaults')

                plugin_image_uri = generate_plugin_image_uri(plugin_name, source, version)
                print(f'Generated plugin_image_uri: {plugin_image_uri}')
            except Exception as e:
                logger.warning(f'Failed to generate plugin_image_uri for {plugin_name}: {e}')
        
        start_task(user_id, task_id, workflow_id, start_time, algorithm_id, plugin_name, task_type, plugin_image_uri, plugin_id)

    def on_success(self, retval, task_id: str, args, kwargs):
        end_time = datetime.now(timezone.utc)
        print(f'Task {task_id} completed at {end_time}, return value: {retval}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='SUCCESS')
        
        # Handle visualization result caching
        task_type = kwargs.get('task_type')
        if task_type == 'visualization':
            cache_key = kwargs.get('cache_key')
            cache_info = kwargs.get('cache_info')

            if cache_key and cache_info:
                try:
                    success = save_result_to_cache(
                        result_file_path=cache_info['result_file_path'],
                        user_path=cache_info['user_path'],
                        cache_key=cache_key,
                        plugin_name=cache_info['plugin_name'],
                        script_name=cache_info['script_name'],
                        linked_location=cache_info['linked_location']
                    )
                    if success:
                        logger.info(f"Successfully cached visualization result for task {task_id}, cache_key: {cache_key}")
                    else:
                        logger.warning(f"Failed to cache visualization result for task {task_id}")
                except Exception as e:
                    logger.error(f"Error caching visualization result for task {task_id}: {e}")
            else:
                logger.debug(f"No cache information provided for visualization task {task_id}")

        if task_type == 'plugin_build':
            try:
                from app.common.utils.plugin_cache import invalidate_all_plugin_cache
                invalidate_all_plugin_cache()
            except Exception as e:
                logger.warning(f"Plugin cache invalidation failed for task {task_id}: {e}")

        # 작업 완료 시 컨테이너 매니저에서 등록 해제
        container_manager.unregister_container(task_id)

    def on_failure(self, exc, task_id: str, args, kwargs, einfo):
        """Ensure the failure is logged and state is correctly updated."""
        logger.error(f"Task {task_id} failed due to {exc}")
        end_time = datetime.now(timezone.utc)
        print(f'Task {task_id} failed at {end_time}, error: {exc}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='FAILURE')

        task_type = kwargs.get('task_type')
        if task_type == 'plugin_build':
            try:
                from app.common.utils.plugin_cache import invalidate_all_plugin_cache
                invalidate_all_plugin_cache()
            except Exception as e:
                logger.warning(f"Plugin cache invalidation failed for task {task_id}: {e}")

        # 작업 실패 시 관련 컨테이너 정리
        try:
            if container_manager.stop_task_container(task_id):
                print(f"Container for failed task {task_id} cleaned up successfully")
            else:
                print(f"No container found or cleanup failed for task {task_id}")
        except Exception as cleanup_error:
            logger.error(f"Error cleaning up container for failed task {task_id}: {cleanup_error}")

        # 작업 실패 시 results 폴더 정리 (재실행 시 Snakemake 스킵 방지)
        try:
            workflow_id = kwargs.get('workflow_id')
            algorithm_id = kwargs.get('algorithm_id')
            task_type = kwargs.get('task_type')

            if user_id and workflow_id and algorithm_id:
                from app.database.conn import get_db_session
                from app.database import models
                from pathlib import Path
                from app.common.utils.log_archive_utils import cleanup_task_results

                with get_db_session() as db:
                    user = db.query(models.User).filter_by(id=user_id).first()
                    if user:
                        if task_type == 'visualization':
                            task_folder = f"./user/{user.username}/workflow_{workflow_id}/visualization_{algorithm_id}"
                        else:
                            task_folder = f"./user/{user.username}/workflow_{workflow_id}/algorithm_{algorithm_id}"

                        cleanup_result = cleanup_task_results(Path(task_folder), preserve_folder=True)
                        if cleanup_result["success"]:
                            files_count = len(cleanup_result.get('files_removed', [])) + len(cleanup_result.get('symlinks_removed', []))
                            print(f"Results cleanup for failed task {task_id}: {files_count} items removed")
                        else:
                            print(f"Results cleanup warning for failed task {task_id}: {cleanup_result.get('error')}")
        except Exception as results_cleanup_error:
            logger.warning(f"Results cleanup failed for task {task_id}: {results_cleanup_error}")

    def on_revoke(self, task_id: str, kwargs, terminated, signum, expired):
        end_time = datetime.now(timezone.utc)
        print(f'Task {task_id} revoked at {end_time}')
        print(f'Revoke details - terminated: {terminated}, signal: {signum}, expired: {expired}')

        user_id = kwargs.get('user_id') if kwargs else None
        if user_id:
            end_task(user_id, task_id, end_time, status='REVOKED')

        task_type = kwargs.get('task_type') if kwargs else None
        if task_type == 'plugin_build':
            try:
                from app.common.utils.plugin_cache import invalidate_all_plugin_cache
                invalidate_all_plugin_cache()
            except Exception as e:
                logger.warning(f"Plugin cache invalidation failed for task {task_id}: {e}")

        # 작업 취소 시 관련 컨테이너 강제 정리
        try:
            print(f"Attempting to stop container for revoked task {task_id}")
            
            # 1. 컨테이너 매니저를 통한 정리
            container_stopped = container_manager.stop_task_container(task_id, timeout=5)
            
            # 2. 컨테이너 이름 패턴으로도 정리 시도 (백업 방법)
            if not container_stopped:
                # task_id 기반 패턴 매칭 (더 정확한 패턴)
                task_short_id = task_id[:8]
                container_patterns = [
                    f"*task-{task_short_id}*",  # 기존 패턴
                    f"plugin-*-task-{task_short_id}-*",  # 더 구체적인 패턴
                ]
                
                for pattern in container_patterns:
                    print(f"Trying to stop containers with pattern: {pattern}")
                    if container_manager.stop_container_by_name(pattern):
                        container_stopped = True
                        break
            
            print(f"Container cleanup completed for revoked task {task_id}")

        except Exception as cleanup_error:
            logger.error(f"Error cleaning up container for revoked task {task_id}: {cleanup_error}")
            print(f"Container cleanup failed for task {task_id}: {cleanup_error}")

        # Snakemake lock 파일 정리 (작업 취소 시 lock이 남아있을 수 있음)
        try:
            snakemake_lock_dir = Path("./").resolve() / ".snakemake" / "locks"
            if snakemake_lock_dir.exists():
                lock_files_removed = 0
                for lock_file in snakemake_lock_dir.glob("*.lock"):
                    lock_file.unlink()
                    lock_files_removed += 1
                if lock_files_removed > 0:
                    print(f"Removed {lock_files_removed} Snakemake lock file(s) for revoked task {task_id}")
        except Exception as lock_cleanup_error:
            logger.warning(f"Snakemake lock cleanup failed for task {task_id}: {lock_cleanup_error}")

    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        print('----------------------------------------')
        print(f'Task {task_id} returned with status {status}, return value: {retval}')
        print('----------------------------------------')
        
        # 모든 작업 완료 후 컨테이너 정리 확인
        if status in ['SUCCESS', 'FAILURE', 'REVOKED']:
            try:
                container_manager.unregister_container(task_id)
            except Exception as e:
                logger.warning(f"Error unregistering container for task {task_id}: {e}")

def wait_for_file_ready(file_path: str, timeout: int = 10, check_interval: float = 0.1) -> bool:
    """
    Wait for a file to be completely written and synced from Docker container.

    Args:
        file_path: Path to the target file
        timeout: Maximum time to wait in seconds (default: 10)
        check_interval: Time between checks in seconds (default: 0.1)

    Returns:
        True if file is ready and stable, False if timeout occurred
    """
    start_time = time.time()
    last_size = -1
    stable_count = 0

    logger.info(f"Waiting for file to be ready: {file_path}")

    while time.time() - start_time < timeout:
        if not Path(file_path).exists():
            time.sleep(check_interval)
            continue

        try:
            current_size = Path(file_path).stat().st_size

            # File size is stable when it remains the same for 3 consecutive checks
            if current_size == last_size and current_size > 0:
                stable_count += 1
                if stable_count >= 3:
                    elapsed = time.time() - start_time
                    logger.info(f"File ready after {elapsed:.2f}s: {file_path} ({current_size} bytes)")
                    return True
            else:
                stable_count = 0

            last_size = current_size
            time.sleep(check_interval)
        except OSError as e:
            # Handle race condition where file exists but stat fails
            logger.debug(f"Temporary error checking file: {e}")
            time.sleep(check_interval)

    # Timeout occurred, but check one last time if file exists
    file_exists = Path(file_path).exists()
    elapsed = time.time() - start_time

    if file_exists:
        logger.warning(f"File exists but may not be stable after {elapsed:.2f}s: {file_path}")
    else:
        logger.error(f"File not found after {elapsed:.2f}s: {file_path}")

    return file_exists


@shared_task(bind=True, base=MyTask, name="workflow_task:process_data_task")
def process_data_task(self, username: str, snakefile_path: str, selected_plugin: str,
                      targets: list, user_id: int, workflow_id: int, algorithm_id: int, plugin_name: str, task_type: str, **kwargs):
    try:
        task_id = self.request.id
        print(f'Processing data for user {username}...')
        print(f"Task ID: {task_id}")
        print(f"Targets: {targets}")
        print(f"Snakefile path: {snakefile_path}")
        print(f"Plugin name: {selected_plugin}")

        self.update_state(state="RUNNING", meta={"message": "Executing workflow..."})

        # Docker 컨테이너로 Snakemake 실행 (task_id 전달)
        result = snakemakeProcess(targets, snakefile_path, selected_plugin, task_id)

        # 실행 결과 검증
        if result["returncode"] != 0:
            error_message = result.get("stderr", "Unknown error occurred")
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)

        # 타겟 파일 존재 여부 확인 (Docker volume sync 대기)
        missing_targets = []
        for target in targets:
            if not wait_for_file_ready(target, timeout=10):
                missing_targets.append(target)

        if missing_targets:
            error_message = f"Target(s) not produced or not synced from container: {missing_targets}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)

        print('Data processing complete.')
        return {
            "status": "Success", 
            "message": "Processing complete",
            "stdout": result.get("stdout", ""),
            "stderr": result.get("stderr", ""),
            "log_path": result.get("log_path", ""),
            "container_id": result.get("container_id", ""),
            "container_name": result.get("container_name", "")
        }

    except Exception as e:
        error_message = str(e)
        if "Plugin image" in error_message:
            error_message = f"Plugin execution failed: {error_message}. Please ensure the plugin is properly built and available."
        self.update_state(state="FAILURE", meta={"error": error_message})
        raise RuntimeError(error_message) from e

@shared_task(bind=True, base=MyTask, name="plugin_task:build_plugin_task")
def build_plugin_task(self, plugin_name: str = None, user_id: int = None, workflow_id: int = None, algorithm_id: int = None, task_type: str = "plugin_build"):
    """
    플러그인 Docker 이미지를 비동기적으로 빌드하는 Celery task
    
    Parameters:
        plugin_name (str): 빌드할 플러그인 이름
        user_id (int): 사용자 ID
        workflow_id (int, optional): 워크플로우 ID (기본값: None)
        algorithm_id (int, optional): 알고리즘 ID (기본값: None)
        task_type (str): 태스크 타입 (기본값: "plugin_build")
    
    Returns:
        dict: 빌드 결과 정보
    """
    try:
        # 필수 매개변수 검증
        if not plugin_name:
            raise ValueError("plugin_name is required")
        if not user_id:
            raise ValueError("user_id is required")
            
        task_id = self.request.id
        print(f'Building plugin Docker image for {plugin_name}...')
        print(f"Task ID: {task_id}")
        print(f"User ID: {user_id}")

        self.update_state(state="RUNNING", meta={"message": f"Building Docker image for plugin {plugin_name}..."})

        # Use new plugin path resolution
        from app.common.utils.plugin_utils import get_plugin_path, is_plugin_editable
        
        # Check if plugin is editable (only local plugins can be built)
        if not is_plugin_editable(plugin_name):
            error_message = f"Cannot build official plugin '{plugin_name}'. Only local plugins can be built."
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)
        
        # Get the plugin path (should be local)
        try:
            plugin_folder, source = get_plugin_path(plugin_name, "local")
        except Exception as e:
            error_message = f"Plugin '{plugin_name}' not found: {str(e)}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)
        
        # 플러그인 폴더가 존재하는지 확인
        if not os.path.exists(plugin_folder):
            error_message = f"Plugin folder not found: {plugin_name}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)
        
        # Dockerfile이 존재하는지 확인
        dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
        if not os.path.exists(dockerfile_path):
            error_message = f"Dockerfile not found in plugin folder: {plugin_name}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)

        # 빌드 상태 업데이트
        self.update_state(state="RUNNING", meta={"message": "Starting Docker build process..."})

        # Docker 이미지 빌드
        build_result = plugin_utils.build_plugin_docker_image(
            plugin_path=plugin_folder,
            plugin_name=plugin_name,
        )

        if not build_result['success']:
            error_message = build_result['message']
            print(error_message)
            self.update_state(state="FAILURE", meta={
                "error": error_message,
                "log_file": build_result.get('log_file'),
                "image_tag": build_result.get('image_tag')
            })
            raise RuntimeError(error_message)

        # Record plugin_image_uri for local build
        try:
            from app.common.utils.plugin_utils import generate_plugin_image_uri
            plugin_image_uri = generate_plugin_image_uri(plugin_name, "local")
            record_plugin_image_uri(task_id, plugin_image_uri, user_id)
            print(f'Recorded plugin_image_uri: {plugin_image_uri}')
        except Exception as e:
            logger.warning(f'Failed to record plugin_image_uri for {plugin_name}: {e}')
        
        print(f'Plugin Docker image build complete for {plugin_name}')
        return {
            "status": "Success", 
            "message": f"Plugin Docker image built successfully for {plugin_name}",
            "plugin_name": plugin_name,
            "log_file": build_result['log_file'],
            "image_tag": build_result['image_tag']
        }

    except Exception as e:
        error_message = str(e)
        if "Docker" in error_message:
            error_message = f"Docker build failed: {error_message}. Please check Docker daemon and plugin configuration."
        self.update_state(state="FAILURE", meta={"error": error_message})
        raise RuntimeError(error_message) from e

