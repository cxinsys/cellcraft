from celery import shared_task, Task
from typing import List
from datetime import datetime
import time
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

    def __call__(self, *args, **kwargs):
        start_time = datetime.now()
        print(f'Task {self.request.id} started at {start_time}')
        user_id = kwargs.get('user_id')
        start_task(user_id, self.request.id, start_time)
        super(MyTask, self).__call__(*args, **kwargs)

def snakemakeProcess(filepath):
    from subprocess import Popen, PIPE
    print(filepath)
    process = Popen(['snakemake',f'workflow/data/{filepath}.csv','-j'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

@shared_task(bind=True, base=MyTask, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5}, name="workflow_task:process_data_task")
def process_data_task(self, username: str, linked_nodes: List[dict], user_id: int):
    # from multiprocessing import Pool, cpu_count
    from billiard import Pool, cpu_count
    print(f'Processing data for user {username}...')
    for nodes in linked_nodes:
        fileName = nodes['file'].replace('.h5ad', '')
        lastNode = "file"
        target = f'{lastNode}_{username}_{fileName}'
        p = Pool(cpu_count())
        snakemake = p.apply_async(snakemakeProcess, (target,))
        print(snakemake.get())
        p.close()
        p.join()
    time.sleep(10)
    print('Data processing complete.')
    return {"status": "Processing complete"}