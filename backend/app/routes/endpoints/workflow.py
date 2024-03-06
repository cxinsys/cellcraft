from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse, StreamingResponse
from sse_starlette.sse import EventSourceResponse
from typing import Any
from sqlalchemy.orm import Session
from celery.result import AsyncResult
from celery.states import REVOKED
import os
import json
import asyncio
from datetime import datetime

from app.common.celery_utils import get_task_info
from app.common.snakemake_utils import create_snakefile, filter_and_add_suffix
from app.database.crud import crud_workflow, crud_task
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
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        user_path = f"./user/{current_user.username}/"
        notAlgoCount = 0
        if user_workflow:
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
            linked_nodes = workflow.linked_nodes
            num = 0
            for nodes in linked_nodes:
                # nodes["lastNode"]가 Algorithm이 아니면 for문을 돌지 않음
                if nodes["lastNode"] != "Algorithm":
                    notAlgoCount += 1
                    continue
                num += 1
                workflow_path = f"{user_path}workflow{num}"
                # user 폴더에 workflow 폴더 존재하지 않으면 생성
                if not os.path.exists(workflow_path):
                    os.makedirs(workflow_path)

                # print(nodes)

                # workflow.option_file(option.json 파일 경로) 가져와서 해당 경로에 있는 json 파일 load
                with open(user_path + "data/" + nodes["algorithmOptions"]["optionFilePath"], "r") as f:
                    options = json.load(f)

                # user_input 데이터를 생성
                user_input = {
                    'userName': current_user.username,  # 현재 사용자 이름
                    'fileName': filter_and_add_suffix(nodes['file']), # 파일 이름
                    'optionFileName': nodes["algorithmOptions"]["optionFilePath"],  # 옵션 파일 이름
                    'numOfThreads': str(options['num_of_threads']), # 스레드 개수
                    'historyLength': str(options['history_length']), # 히스토리 길이
                    'cutoffForFdr': str(options['cutoff_for_fdr']),  # FDR 기준
                    'numOfLinks': str(options['num_of_links']), # 링크 개수
                    'species': options['species'], # 종
                    'trimmingIndirectEdges': str(options['trimming_indirect_edges'])  # 간접 에지 자르기
                }
                # user 폴더에 snakefile 생성
                snakefile_path = create_snakefile(workflow_path + "/snakefile", user_input)

                # target을 nodes.algorithmOptions.algorithm이 TENET일 때랑 TENET_TF일 때로 나눠서 처리
                if nodes["algorithmOptions"]["algorithm"] == "TENET":
                    target = f"workflow/data/DownstreamAnalysis_{current_user.username}_{filter_and_add_suffix(nodes['file'])}.txt"
                elif nodes["algorithmOptions"]["algorithm"] == "TENET_TF":
                    target = f"workflow/data/DownstreamAnalysisTF_{current_user.username}_{filter_and_add_suffix(nodes['file'])}.txt"

                process_task = process_data_task.apply_async(
                    (current_user.username, snakefile_path, target),
                    kwargs={'user_id': current_user.id, 'workflow_id': workflow.id }
                )
                message = "Tasks added to queue"
                task_id = process_task.id
                result = get_task_info(process_task.id)
        else:
            workflow_info = crud_workflow.create_workflow(db, workflow.title, workflow.thumbnail, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)
            linked_nodes = workflow.linked_nodes
            num = 0
            for nodes in linked_nodes:
                # nodes["lastNode"]가 Algorithm이 아니면 for문을 돌지 않음
                if nodes["lastNode"] != "Algorithm":
                    notAlgoCount += 1
                    continue
                num += 1
                workflow_path = f"{user_path}workflow{num}"
                # user 폴더에 workflow 폴더 존재하지 않으면 생성
                if not os.path.exists(workflow_path):
                    os.makedirs(workflow_path)

                # print(nodes)

                # workflow.option_file(option.json 파일 경로) 가져와서 해당 경로에 있는 json 파일 load
                with open(user_path + "data/" + nodes["algorithmOptions"]["optionFilePath"], "r") as f:
                    options = json.load(f)

                # user_input 데이터를 생성
                user_input = {
                    'userName': current_user.username,  # 현재 사용자 이름
                    'fileName': filter_and_add_suffix(nodes['file']), # 파일 이름
                    'optionFileName': nodes["algorithmOptions"]["optionFilePath"],  # 옵션 파일 이름
                    'numOfThreads': str(options['num_of_threads']), # 스레드 개수
                    'historyLength': str(options['history_length']), # 히스토리 길이
                    'cutoffForFdr': str(options['cutoff_for_fdr']),  # FDR 기준
                    'numOfLinks': str(options['num_of_links']), # 링크 개수
                    'species': options['species'], # 종
                    'trimmingIndirectEdges': str(options['trimming_indirect_edges'])  # 간접 에지 자르기
                }
                # user 폴더에 snakefile 생성
                snakefile_path = create_snakefile(workflow_path + "/snakefile", user_input)

                # target을 nodes.algorithmOptions.algorithm이 TENET일 때랑 TENET_TF일 때로 나눠서 처리
                if nodes["algorithmOptions"]["algorithm"] == "TENET":
                    target = f"workflow/data/DownstreamAnalysis_{current_user.username}_{filter_and_add_suffix(nodes['file'])}.txt"
                elif nodes["algorithmOptions"]["algorithm"] == "TENET_TF":
                    target = f"workflow/data/DownstreamAnalysisTF_{current_user.username}_{filter_and_add_suffix(nodes['file'])}.txt"

                process_task = process_data_task.apply_async(
                    (current_user.username, snakefile_path, target),
                    kwargs={'user_id': current_user.id, 'workflow_id': workflow.id }
                )
                message = "Tasks added to queue"
                task_id = process_task.id
                result = get_task_info(process_task.id)
        
        if notAlgoCount == len(linked_nodes):
            # 에러 일으키기
            raise HTTPException(
                status_code=400,
                detail="please select last node as Algorithm",
                )

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
        return crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
    else :
        # workflow 생성
        return crud_workflow.create_workflow(db, workflow.title, workflow.thumbnail, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)

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
                'thumbnail': item.thumbnail, # 시현 추가
                'updated_at': item.updated_at, 
                'user_id': item.user_id
            }
            res.append(workflow_res)
            # print(res)
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
            'thumbnail': user_workflow.thumbnail, # 시현 추가
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
    # print(file_list)
    FILE_NAME = filename.filename
    
    for item_file in file_list:
        if FILE_NAME in item_file:
            FILE_NAME = item_file
    # print(FILE_NAME)
    FILE_PATH = os.path.join(PATH_COMPILE_RESULT, FILE_NAME)
    # print(FILE_PATH)

    return FileResponse(FILE_PATH)

#response workflow result
@router.post("/results")
def getResults(current_user: models.User = Depends(dep.get_current_active_user)):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/result'
    file_list = os.listdir(PATH_COMPILE_RESULT)

    return file_list

@router.get("/task/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    async def event_generator():
        while True:
            if task_id:
                task = get_task_info(task_id)
                if task['task_status'] == 'SUCCESS' or task['task_status'] == 'FAILURE' or task['task_status'] == 'REVOKED' or task['task_status'] == 'RETRY':
                    yield f"{task['task_status']}"
                    break
                print(task['task_status'])
                yield f"{task['task_status']}"
                await asyncio.sleep(5)
            else:
                break
    return EventSourceResponse(event_generator())

@router.get("/monitoring")
async def get_task_monitoring(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any :
    """
    Return the status of the all User Task
    """
    user_task = crud_task.get_user_task(db, current_user.id)
    if user_task:
        return user_task
    else:
        raise HTTPException(
                status_code=400,
                detail="this user not exists task",
                )

@router.delete("/revoke/{task_id}")
def revoke_task(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
    ) -> dict:
    """
    Revoke the task
    """
    from app.main import get_celery_app
    celery = get_celery_app()
    celery.control.revoke(task_id, terminate=True, signal='SIGTERM')
    task = get_task_info(task_id)

    print(task)

    # task_info가 사전 형식인 경우 상태를 접근하는 방식
    task_status = task.get("status")
    if task_status == 'REVOKED':
        return {"message": "Task Revoked", "task_id": task_id}
    else:
        # 태스크 상태를 'REVOKED'로 업데이트
        crud_task.end_task(current_user.id, task_id, datetime.now(), 'REVOKED')
        # 태스크가 업데이트 되었는지 확인
        task = get_task_info(task_id)
        task_status = task.get("status")
        if task_status == 'REVOKED':
            return {"message": "Task Revoked", "task_id": task_id}
        else:
            return {"message": "Task Revoked Failed", "task_id": task_id}
