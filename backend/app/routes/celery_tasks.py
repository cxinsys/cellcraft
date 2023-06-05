from celery import shared_task
from typing import List
from app.database.schemas.workflow import WorkflowCreate

def snakemakeProcess(filepath):
    from subprocess import Popen, PIPE
    print(filepath)
    process = Popen(['snakemake',f'workflow/data/{filepath}.csv','-j'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5}, name="workflow_task:process_data_task")
def process_data_task(self, username: str, linked_nodes: List[dict]):
    # from multiprocessing import Pool, cpu_count
    from billiard import Pool, cpu_count
    for nodes in linked_nodes:
        fileName = nodes['file'].replace('.h5ad', '')
        lastNode = "file"
        target = f'{lastNode}_{username}_{fileName}'
        p = Pool(cpu_count())
        snakemake = p.apply_async(snakemakeProcess, (target,))
        print(snakemake.get())
        p.close()
        p.join()
    return {"status": "Processing complete"}

@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5}, name="workflow_db_task:manage_workflow_task")
def manage_workflow_task(self, db, user_id: int, workflow: WorkflowCreate):
    from . import crud_workflow  # import here to avoid circular imports
    user_workflow = crud_workflow.get_user_workflow(db, user_id, workflow.id)
    if user_workflow:
        crud_workflow.update_workflow(db, user_id, workflow.id, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
    else:
        crud_workflow.create_workflow(db, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, user_id)
    return {"status": "Workflow management complete"}
