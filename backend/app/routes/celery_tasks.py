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

# 통합된 리소스 카운터
class ResourceCounter:
    def __init__(self, max_total_tasks, max_gpu_tasks):
        self.max_total_tasks = max_total_tasks
        self.max_gpu_tasks = max_gpu_tasks
        self.current_total_tasks = 0
        self.current_gpu_tasks = 0
        self._lock = threading.Lock()
    
    def acquire(self, use_gpu=False):
        with self._lock:
            if self.current_total_tasks >= self.max_total_tasks:
                return False
                
            if use_gpu:
                if self.current_gpu_tasks >= self.max_gpu_tasks:
                    return False
                self.current_gpu_tasks += 1
            
            self.current_total_tasks += 1
            return True
    
    def release(self, use_gpu=False):
        with self._lock:
            if use_gpu:
                self.current_gpu_tasks = max(0, self.current_gpu_tasks - 1)
            self.current_total_tasks = max(0, self.current_total_tasks - 1)

# 통합된 리소스 카운터 인스턴스 생성
resource_counter = ResourceCounter(max_total_tasks=11, max_gpu_tasks=7)

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
    
    cpu_cores = 1 if use_gpu else 4  # GPU 태스크는 1코어, CPU 태스크는 4코어 사용

    retry_count = 0
    while retry_count < 50:
        if resource_counter.acquire(use_gpu):
            try:
                print(f'Processing data for user {username}...')
                print(f'Using {cpu_cores} CPU cores, GPU: {use_gpu}')

                # 의존성 설치
                for dependency_file in ['requirements.txt', 'environment.yml', 'environment.yaml', 'renv.lock']:
                    dependency_file_path = os.path.join(plugin_dependency_path, dependency_file)
                    if os.path.exists(dependency_file_path):
                        print(f"Installing dependencies from {dependency_file}...")
                        try:
                            install_dependencies(dependency_file_path)
                        except Exception as e:
                            # 의존성 설치 실패 시 태스크 상태를 FAILURE로 설정
                            error_message = f"Failed to install dependencies from {dependency_file}: {str(e)}"
                            print(error_message)
                            raise RuntimeError(error_message)  # 태스크 실패


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
                resource_counter.release(use_gpu)
                break
        else:
            retry_count += 1
            end_time = datetime.now()
            end_task(user_id, self.request.id, end_time, status='PENDING')
            self.retry(countdown=10, max_retries=None)

    if retry_count >= 50:
        raise RuntimeError("Failed to acquire resources after multiple attempts")
