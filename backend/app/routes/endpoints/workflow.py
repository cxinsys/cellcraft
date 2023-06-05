from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from typing import Any
from sqlalchemy.orm import Session
import os
import json

from app.common.celery_utils import get_task_info
from app.database.crud import crud_workflow
from app.database.schemas.workflow import WorkflowDelete, WorkflowCreate, WorkflowResult, WorkflowFind
from app.routes import dep
from app.database import models
from app.routes.celery_tasks import process_data_task

router = APIRouter()                                                        

#export workflow data
@router.post("/compile")
def exportData(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        process_task = process_data_task.apply_async(
            (current_user.username, workflow.linked_nodes),
            kwargs={'user_id': current_user.id}
        )
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        if user_workflow:
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
        else:
            crud_workflow.create_workflow(db, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)
        message = "Tasks added to queue"
        task_id = process_task.id
        result = get_task_info(process_task.id)

    except json.JSONDecodeError:
        workflow = None
        message = "Received data is not a valid JSON"
        task_id = {}
    return {"message": message, "recived_data": workflow, "task_id": task_id, "result": result}


#save workflow data
@router.post("/save")
def update_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
    if user_workflow:
        # workflow 수정
        return crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
    else :
        # workflow 생성
        return crud_workflow.create_workflow(db, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)

#User workflow delete
@router.post("/delete")
def delete_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowDelete,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        delete_workflow = crud_workflow.delete_user_workflow(db, current_user.id, workflowInfo.id)
        print(delete_workflow)
        return delete_workflow
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#response workflow data
@router.get("/me")
def get_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    user_workflows = crud_workflow.get_user_workflows(db, current_user.id)
    if user_workflows:
        res = []
        for item in user_workflows:
            workflow_res = {
                'id': item.id, 
                'title': item.title, 
                'updated_at': item.updated_at, 
                'user_id': item.user_id
            }
            res.append(workflow_res)
            print(res)
        return res
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#find user workflow
@router.post("/find", response_model=WorkflowCreate)
def find_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowFind,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        return {
            'title': user_workflow.title,
            'workflow_info': user_workflow.workflow_info,
            'nodes': user_workflow.nodes,
            'linked_nodes': user_workflow.linked_nodes,
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )


#response workflow result
@router.post("/result")
def checkResult(filename: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/result'
    file_list = os.listdir(PATH_COMPILE_RESULT)
    print(file_list)
    FILE_NAME = filename.filename
    
    for item_file in file_list:
        if FILE_NAME in item_file:
            FILE_NAME = item_file
    print(FILE_NAME)
    # FILE_PATH = PATH_COMPILE_RESULT + '/' + FILE_NAME 
    FILE_PATH = PATH_COMPILE_RESULT + '/' + "file_pbmc3k_obs_umap.csv"

    # if FILE_NAME.find("Plot") != -1:
    #     # img = Image.open(FILE_PATH)
    #     # img_converted = from_image_to_bytes(img)
    #     # return JSONResponse(img_converted)
    #     with open(FILE_PATH, 'r') as file:
    #         data = json.load(file)
    #         print(data)
    #         return JSONResponse(data)
    # else:
    #     return FileResponse(FILE_PATH)
    return FileResponse(FILE_PATH)

@router.get("/task/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    return get_task_info(task_id)