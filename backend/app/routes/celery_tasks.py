from celery import shared_task, Task
from datetime import datetime
from celery.signals import worker_process_init
import os
import time
from typing import List
from billiard import Pool, cpu_count
import amqp

from app.common.utils.snakemake_utils import snakemakeProcess
from app.common.utils.plugin_utils import install_dependencies
from app.database.crud.crud_task import start_task, end_task

class MyTask(Task):
    def on_success(self, retval, task_id: str, args, kwargs):
        end_time = datetime.now()
        print(f'Task {task_id} completed at {end_time}, return value: {retval}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='SUCCESS')

    def on_failure(self, exc, task_id: str, args, kwargs, einfo):
        end_time = datetime.now()
        print(f'Task {task_id} failed at {end_time}, error: {exc}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='FAILURE')

    def on_revoke(self, task_id: str, kwargs, terminated, signum, expired):
        end_time = datetime.now()
        print(f'Task {task_id} revoked at {end_time}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='REVOKED')

    def __call__(self, *args, **kwargs):
        start_time = datetime.now()
        print(f'Task {self.request.id} started at {start_time}')
        user_id = kwargs.get('user_id')
        workflow_id = kwargs.get('workflow_id')
        start_task(user_id, self.request.id, workflow_id, start_time)
        super(MyTask, self).__call__(*args, **kwargs)

@worker_process_init.connect
def configure_worker(conf=None, **kwargs):
    os.environ['CONDA_DEFAULT_ENV'] = 'snakemake'

@shared_task(bind=True, base=MyTask, name="workflow_task:process_data_task")
def process_data_task(self, username: str, snakefile_path: str, plugin_dependency_path: str, targets: List[str], user_id: int, workflow_id: int):

    try:
        print(f'Processing data for user {username}...')
        print(f'Creating a new process pool with {cpu_count()} processes...')
        print(f'Targets: {targets}')
        print(f'Snakefile path: {snakefile_path}')

        # 의존성 설치
        for dependency_file in ['requirements.txt', 'environment.yml', 'environment.yaml', 'renv.lock']:
            dependency_file_path = os.path.join(plugin_dependency_path, dependency_file)
            if os.path.exists(dependency_file_path):
                print(f"Installing dependencies from {dependency_file}...")
                install_dependencies(dependency_file_path)

        p = Pool(cpu_count())
        snakemake_process = p.apply_async(snakemakeProcess, (targets, snakefile_path))
        
        while not snakemake_process.ready():
            try:
                # Celery 작업의 취소 상태 확인
                revoked_tasks = self.app.control.inspect().revoked()
                if revoked_tasks is not None and self.request.id in revoked_tasks.get(self.request.id, []):
                    print(f"Task {self.request.id} was revoked.")
                    snakemake_process.terminate()
                    break
            except Exception as e:
                print(f"오류 발생: {e}")

        result = snakemake_process.get()
        p.close()
        p.join()

        if result["returncode"] != 0:
            raise RuntimeError(f"Snakemake process failed: {result['stderr']}")

        print('Data processing complete.')
        return {"status": "Success", "message": "Processing complete"}

    except amqp.exceptions.PreconditionFailed as e:
        print(f"메시지 브로커 연결 오류 발생: {e}")
        # 메시지 브로커에 재연결 시도
        try:
            self.app.connection().ensure_connection(max_retries=3)
            print("메시지 브로커에 재연결 시도 성공")
        except Exception as e:
            print(f"재연결 시도 실패: {e}")
            raise e
    except Exception as e:
        # 오류 발생 시 on_failure 메서드가 자동으로 호출됨
        raise e