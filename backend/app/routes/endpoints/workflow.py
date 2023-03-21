from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File, responses
from fastapi.responses import FileResponse, JSONResponse
from fastapi.encoders import jsonable_encoder
from typing import List, Union, Any
from subprocess import Popen,PIPE
from multiprocessing import Process, Pool, cpu_count
from sqlalchemy.orm import Session
import os
import json
import base64
import io
from PIL import Image

from app.database.crud import crud_workflow
from app.database.schemas.workflow import WorkflowCreate, WorkflowResult, WorkflowFind
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
        crud_workflow.create_workflow(db, workflow.title, workflow.workflow_info, workflow.nodes, workflow.linked_nodes, current_user.id)
        message = "success"
    except json.JSONDecodeError:
        workflow = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": workflow}

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