from celery import shared_task, Task
from datetime import datetime
from celery.signals import worker_process_init
import os
import time
from typing import List
from billiard import Pool, cpu_count
import amqp
from functools import partial
import threading

from app.common.utils.snakemake_utils import snakemakeProcess
from app.common.utils.plugin_utils import install_dependencies
from app.database.crud.crud_task import start_task, end_task

# 전역 리소스 카운터
class ResourceCounter:
    def __init__(self, max_tasks):
        self.max_tasks = max_tasks
        self.current_tasks = 0
        self._lock = threading.Lock()
    
    def acquire(self):
        with self._lock:
            if self.current_tasks < self.max_tasks:
                self.current_tasks += 1
                return True
            return False
    
    def release(self):
        with self._lock:
            self.current_tasks = max(0, self.current_tasks - 1)

# CPU 워커용 카운터 (44코어 / 4코어 = 최대 11개 태스크)
cpu_resource_counter = ResourceCounter(max_tasks=11)
# GPU 워커용 카운터 (7개 GPU = 최대 7개 태스크)
gpu_resource_counter = ResourceCounter(max_tasks=7)

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
def process_data_task(self, username: str, snakefile_path: str, plugin_dependency_path: str, 
                     targets: List[str], user_id: int, workflow_id: int, 
                     use_gpu: bool = False):
    
    # GPU/CPU 태스크 구분
    if use_gpu:
        resource_counter = gpu_resource_counter
        cpu_cores = 1  # GPU 태스크는 1개 CPU 코어 사용
        gpu_cores = 1
        self.request.delivery_info['routing_key'] = 'gpu_tasks'
    else:
        resource_counter = cpu_resource_counter
        cpu_cores = 4  # CPU 태스크는 4개 코어 사용
        gpu_cores = 0
        self.request.delivery_info['routing_key'] = 'cpu_tasks'

    retry_count = 0
    while retry_count < 3:
        if resource_counter.acquire():
            try:
                print(f'Processing data for user {username}...')
                print(f'Using {cpu_cores} CPU cores and {gpu_cores} GPU cores')

                # 의존성 설치
                for dependency_file in ['requirements.txt', 'environment.yml', 'environment.yaml', 'renv.lock']:
                    dependency_file_path = os.path.join(plugin_dependency_path, dependency_file)
                    if os.path.exists(dependency_file_path):
                        print(f"Installing dependencies from {dependency_file}...")
                        install_dependencies(dependency_file_path)

                # CPU 코어 수에 맞게 Pool 생성
                p = Pool(cpu_cores)
                snakemake_process = p.apply_async(snakemakeProcess, (targets, snakefile_path))
                
                while not snakemake_process.ready():
                    try:
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

            finally:
                resource_counter.release()
                break
        else:
            retry_count += 1
            end_time = datetime.now()
            end_task(user_id, self.request.id, end_time, status='PENDING')
            self.retry(countdown=60, max_retries=None)

    if retry_count >= 3:
        raise RuntimeError("Failed to acquire resources after multiple attempts")
