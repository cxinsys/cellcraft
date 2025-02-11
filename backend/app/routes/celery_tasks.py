from celery import shared_task, Task
from datetime import datetime
from celery.signals import worker_process_init
from celery.exceptions import Ignore
from celery.worker.request import Request
import os
import time
from typing import List
from billiard import Pool, cpu_count
import amqp
from functools import partial
import threading
import logging
from celery.exceptions import MaxRetriesExceededError

from app.common.utils.snakemake_utils import snakemakeProcess
from app.common.utils.plugin_utils import install_dependencies
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

@worker_process_init.connect
def configure_worker(conf=None, **kwargs):
    os.environ['CONDA_DEFAULT_ENV'] = 'snakemake'

@shared_task(bind=True, base=MyTask, name="workflow_task:process_data_task")
def process_data_task(self, username: str, snakefile_path: str, plugin_dependency_path: str, 
                      targets: list, user_id: int, workflow_id: int):

    try:
        print(f'Processing data for user {username}...')
        print(f"Task ID: {self.request.id}")

        self.update_state(state="INSTALLING", meta={"message": "Installing dependencies..."})

        # 의존성 설치
        for dependency_file in ['requirements.txt', 'environment.yml', 'environment.yaml', 'renv.lock']:
            dependency_file_path = os.path.join(plugin_dependency_path, dependency_file)
            if os.path.exists(dependency_file_path):
                print(f"Installing dependencies from {dependency_file}...")
                try:
                    install_dependencies(dependency_file_path)
                except Exception as e:
                    error_message = f"Failed to install dependencies from {dependency_file}: {str(e)}"
                    print(error_message)
                    self.update_state(state="FAILURE", meta={"error": error_message})
                    raise RuntimeError(error_message)

        self.update_state(state="RUNNING", meta={"message": "Executing workflow..."})

        # Snakemake 실행
        result = snakemakeProcess(targets, snakefile_path)

        print(f"Snakemake process return code: {result['returncode']}")

        # Snakemake 실패 감지 후 Celery 태스크 실패 처리
        if result["returncode"] != 0:
            error_message = f"Snakemake process failed with error:\n{result.get('stderr', 'No error message')}"
            print(error_message)
            self.update_state(state="FAILURE", meta={"error": error_message})

            print("Debug: Raising RuntimeError...")
            raise RuntimeError(error_message)
            
        else:
            print('Data processing complete.')
            retval = {"status": "Success", "message": "Processing complete"}
            return retval

    except Exception as e:
        self.update_state(state="FAILURE", meta={"error": str(e)})
        raise RuntimeError(str(e)) from e

