from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File, responses
from fastapi.responses import FileResponse, JSONResponse
from fastapi.encoders import jsonable_encoder
from typing import List, Union
from subprocess import Popen,PIPE
from multiprocessing import Process, Pool, cpu_count
import os
import json
import base64
import requests
import io
from PIL import Image

from app.database.schemas.workflow import WorkflowResult
from app.routes import dep
from app.database import models

router = APIRouter()

def snakemakeProcess(filepath):
    print(filepath)
    process = Popen(['snakemake',f'workflow/data/{filepath}.csv','-j'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

#export workflow data
@router.post("/compile")
async def exportData(request: Request, current_user: models.User = Depends(dep.get_current_active_user)):
    try:
        payload_as_json = await request.json()
        # print(type(payload_as_json))
        for nodes in payload_as_json:
            # print(type(nodes))
            fileName = nodes['file'].replace('.csv', '')
            lastNode = nodes['lastNode']
            print(fileName, lastNode)
            target = f'{lastNode}_{current_user.username}_{fileName}'
            p = Pool(cpu_count())
            snakemake = p.apply_async(snakemakeProcess, (target,))
            print(snakemake.get())
            p.close()
            p.join()
        message = "success"
    except json.JSONDecodeError:
        payload_as_json = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": payload_as_json}

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
    FILE_PATH = PATH_COMPILE_RESULT + '/' + FILE_NAME 

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