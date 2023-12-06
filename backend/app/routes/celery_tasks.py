from celery import shared_task, Task
from typing import List
from datetime import datetime
from celery.signals import worker_process_init
import time
import os
import signal
import json


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

    def on_retry(self, exc, task_id: str, args, kwargs, einfo):
        end_time = datetime.now()
        print(f'Task {task_id} retried at {end_time}, error: {exc}')
        user_id = kwargs.get('user_id')
        end_task(user_id, task_id, end_time, status='RETRY')

    def __call__(self, *args, **kwargs):
        start_time = datetime.now()
        print(f'Task {self.request.id} started at {start_time}')
        user_id = kwargs.get('user_id')
        workflow_id = kwargs.get('workflow_id')
        start_task(user_id, self.request.id, workflow_id, start_time)
        super(MyTask, self).__call__(*args, **kwargs)

def snakemakeProcess(filepath):
    from subprocess import Popen, PIPE
    print(filepath)
    print(os.environ["PATH"])

    process = Popen(['/opt/conda/envs/snakemake/bin/snakemake', f'workflow/data/{filepath}.txt', '-j'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)


def is_snakemake_installed() -> bool:
    import subprocess
    try:
        # Try running `snakemake --version`
        subprocess.run(["snakemake", "--version"], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return False
    
def list_conda_libraries():
    import subprocess
    try:
        result = subprocess.run(["conda", "list"], check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching installed libraries: {e}")
        return None
    
def list_pip_libraries():
    import subprocess
    try:
        result = subprocess.run(["pip", "list"], check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching installed libraries: {e}")
        return None

def get_conda_environment_name():
    import sys
    print(sys.executable)
    stat_info = os.stat('/app')
    owner_uid = stat_info.st_uid
    owner_gid = stat_info.st_gid

    print(f"Owner UID: {owner_uid}")
    print(f"Owner GID: {owner_gid}")

    return os.environ.get('CONDA_DEFAULT_ENV', None)

def filter_and_add_suffix(input_string):
    # Check if the input_string contains an underscore
    if "_" in input_string:
        # Find the position of the first underscore
        underscore_position = input_string.index("_")
        # Add ".h5ad" before the first underscore and exclude everything after underscore
        modified_string = input_string[:underscore_position] + ".h5ad"
        # Return the modified string
        return modified_string
    # If no underscore found, return the original string
    return input_string

@worker_process_init.connect
def configure_worker(conf=None, **kwargs):
    os.environ['CONDA_DEFAULT_ENV'] = 'snakemake'

@shared_task(bind=True, base=MyTask, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5}, name="workflow_task:process_data_task")
def process_data_task(self, username: str, linked_nodes: List[dict], user_id: int, workflow_id: int):
    # from multiprocessing import Pool, cpu_count
    from billiard import Pool, cpu_count
    print(f'Processing data for user {username}...')
    for nodes in linked_nodes:
        fileName = filter_and_add_suffix(nodes['file'])
        file_path = f"user/{username}/data/{fileName}_option.json"
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)

        # 'algorithm' 키의 값을 추출
        algorithm = data.get('algorithm', None)
        if algorithm == 'fasttenet':
            lastNode = 'FastTenet'
        elif algorithm == 'tenet':
            lastNode = 'DownstreamAnalysis'
        else:
            print("'algorithm' 키가 존재하지 않습니다.")
        # lastNode = nodes['lastNode']
        target = f'{lastNode}_{username}_{fileName}'
        print(target)
        p = Pool(cpu_count())
        snakemake_process = p.apply_async(snakemakeProcess, (target,))
        process = snakemake_process.get()

        # Monitor for task termination request
        while process.poll() is None:  # While the process is running
            if self.request.called_directly:  # This checks if the task is being called directly, not by a worker.
                break
            if self.request.terminate:  # Check for termination request
                # If termination is requested, send SIGTERM to snakemake process
                process.send_signal(signal.SIGTERM)
                break

        p.close()
        p.join()
    print('Data processing complete.')
    return {"status": "Processing complete"}