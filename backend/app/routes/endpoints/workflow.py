from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from typing import Any
from subprocess import Popen,PIPE
from multiprocessing import Pool, cpu_count
from sqlalchemy.orm import Session
import os
import json
import base64
import io

from app.database.crud import crud_workflow
from app.database.schemas.workflow import WorkflowDelete, WorkflowCreate, WorkflowResult, WorkflowFind
from app.routes import dep
from app.database import models

router = APIRouter()                                                        

def snakemakeProcess(filepath):
    print(filepath)
    process = Popen(['snakemake',f'workflow/data/{filepath}.csv','-j'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

# 이미지 base64 변환
def from_image_to_bytes(img):
    # Pillow 이미지 객체를 Bytes로 변환
    imgByteArr = io.BytesIO()
    img.save(imgByteArr, format=img.format)
    imgByteArr = imgByteArr.getvalue()
    # Base64로 Bytes를 인코딩
    encoded = base64.b64encode(imgByteArr)
    # Base64로 ascii로 디코딩
    decoded = encoded.decode('ascii')
    return decoded

#export workflow data
@router.post("/compile")
async def exportData(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        # print(workflow.workflow_info)
        # print(workflow.nodes)
        print(workflow.linked_nodes)
        for nodes in workflow.linked_nodes:
            # fileName = nodes['file'].replace('.csv', '')
            fileName = nodes['file'].replace('.h5ad', '')
            # lastNode = nodes['lastNode']
            lastNode = "file"
            print(fileName, lastNode)
            target = f'{lastNode}_{current_user.username}_{fileName}'
            p = Pool(cpu_count())
            snakemake = p.apply_async(snakemakeProcess, (target,))
            print(snakemake.get())
            p.close()
            p.join()
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        if user_workflow:
            # workflow 수정
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes)
        else :
            # workflow 생성
            crud_workflow.create_workflow(db, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)
        message = "success"
    except json.JSONDecodeError:
        workflow = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": workflow}

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