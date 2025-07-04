from celery import shared_task, Task
from datetime import datetime
from celery.signals import worker_process_init
from celery.exceptions import Ignore
from celery.worker.request import Request
import os
from pathlib import Path
import time
from typing import List
from billiard import Pool, cpu_count
import amqp
from functools import partial
import threading
import logging
from celery.exceptions import MaxRetriesExceededError

from app.common.utils.snakemake_utils import snakemakeProcess
from app.common.utils.docker_utils import container_manager
from app.database.crud.crud_task import start_task, end_task

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
        start_time = datetime.now()
        print(f'Task {task_id} started at {start_time}')
        user_id = kwargs.get('user_id')
        workflow_id = kwargs.get('workflow_id')
        algorithm_id = kwargs.get('algorithm_id')
        plugin_name = kwargs.get('plugin_name')
        task_type = kwargs.get('task_type')
        start_task(user_id, task_id, workflow_id, start_time, algorithm_id, plugin_name, task_type)

    def on_success(self, retval, task_id: str, args, kwargs):
        end_time = datetime.now()
        print(f'Task {task_id} completed at {end_time}, return value: {retval}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='SUCCESS')
        
        # 작업 완료 시 컨테이너 매니저에서 등록 해제
        container_manager.unregister_container(task_id)

    def on_failure(self, exc, task_id: str, args, kwargs, einfo):
        """Ensure the failure is logged and state is correctly updated."""
        logger.error(f"Task {task_id} failed due to {exc}")
        end_time = datetime.now()
        print(f'Task {task_id} failed at {end_time}, error: {exc}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='FAILURE')
        
        # 작업 실패 시 관련 컨테이너 정리
        try:
            if container_manager.stop_task_container(task_id):
                print(f"Container for failed task {task_id} cleaned up successfully")
            else:
                print(f"No container found or cleanup failed for task {task_id}")
        except Exception as cleanup_error:
            logger.error(f"Error cleaning up container for failed task {task_id}: {cleanup_error}")

    def on_revoke(self, task_id: str, kwargs, terminated, signum, expired):
        end_time = datetime.now()
        print(f'Task {task_id} revoked at {end_time}')
        print(f'Revoke details - terminated: {terminated}, signal: {signum}, expired: {expired}')
        
        user_id = kwargs.get('user_id') if kwargs else None
        if user_id:
            end_task(user_id, task_id, end_time, status='REVOKED')
        
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

@shared_task(bind=True, base=MyTask, name="workflow_task:process_data_task")
def process_data_task(self, username: str, snakefile_path: str, selected_plugin: str, 
                      targets: list, user_id: int, workflow_id: int, algorithm_id: int, plugin_name: str, task_type: str):
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

        # 타겟 파일 존재 여부 확인
        missing_targets = []
        for target in targets:
            if not Path(target).exists():
                missing_targets.append(target)

        if missing_targets:
            error_message = f"Target(s) not produced: {missing_targets}"
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

