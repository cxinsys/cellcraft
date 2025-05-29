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
        start_task(user_id, task_id, workflow_id, start_time)

    def on_success(self, retval, task_id: str, args, kwargs):
        end_time = datetime.now()
        print(f'Task {task_id} completed at {end_time}, return value: {retval}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='SUCCESS')

    def on_failure(self, exc, task_id: str, args, kwargs, einfo):
        """Ensure the failure is logged and state is correctly updated."""
        logger.error(f"Task {task_id} failed due to {exc}")
        end_time = datetime.now()
        print(f'Task {task_id} failed at {end_time}, error: {exc}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='FAILURE')

    def on_revoke(self, task_id: str, kwargs, terminated, signum, expired):
        end_time = datetime.now()
        print(f'Task {task_id} revoked at {end_time}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='REVOKED')

    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        print('----------------------------------------')
        print(f'Task {task_id} returned with status {status}, return value: {retval}')
        print('----------------------------------------')

    # def __call__(self, *args, **kwargs):
    #     start_time = datetime.now()
    #     print(f'Task {self.request.id} started at {start_time}')
    #     user_id = kwargs.get('user_id')
    #     workflow_id = kwargs.get('workflow_id')
    #     start_task(user_id, self.request.id, workflow_id, start_time)
    #     super().__call__(*args, **kwargs)

@shared_task(bind=True, base=MyTask, name="workflow_task:process_data_task")
def process_data_task(self, username: str, snakefile_path: str, plugin_name: str, 
                      targets: list, user_id: int, workflow_id: int):
    try:
        print(f'Processing data for user {username}...')
        print(f"Task ID: {self.request.id}")
        print(f"Targets: {targets}")
        print(f"Snakefile path: {snakefile_path}")

        self.update_state(state="RUNNING", meta={"message": "Executing workflow..."})

        # 작업 디렉토리 설정
        workspace_path = Path(os.path.dirname(snakefile_path))
        log_dir = workspace_path / "logs"
        log_file = log_dir / "run.log"

        # Docker 컨테이너로 Snakemake 실행
        result = snakemakeProcess(targets, snakefile_path, plugin_name)

        # 로그 파일 검증
        if not log_file.exists() or log_file.stat().st_size == 0:
            error_message = "run.log not created — Snakemake may not have run"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)

        # 타겟 파일 존재 여부 확인
        missing_targets = []
        for target in targets:
            if not target.exists():
                missing_targets.append(target)

        if missing_targets:
            error_message = f"Target(s) not produced: {missing_targets}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})
            raise RuntimeError(error_message)

        print('Data processing complete.')
        retval = {
            "status": "Success", 
            "message": "Processing complete",
            "stdout": result.get("stdout", ""),
            "stderr": result.get("stderr", ""),
            "log_path": str(log_file)
        }
        return retval

    except Exception as e:
        error_message = str(e)
        if "Plugin image" in error_message:
            error_message = f"Plugin execution failed: {error_message}. Please ensure the plugin is properly built and available."
        self.update_state(state="FAILURE", meta={"error": error_message})
        raise RuntimeError(error_message) from e

