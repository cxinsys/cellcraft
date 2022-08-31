from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File, responses
from fastapi.responses import FileResponse, JSONResponse
from fastapi.encoders import jsonable_encoder
from typing import List, Union
from subprocess import Popen,PIPE
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

def linkList(nodeList):
    result = []
    #while문으로 검사
    while True:
        for item in nodeList:
            for i in range(len(nodeList)):
                if(i >= len(nodeList)):
                    continue
                if item[len(item)-1] == nodeList[i][0]:
                    nodeList.append(item + nodeList[i])
                    del nodeList[i]
                    nodeList.remove(item)
        break
    for itemList in nodeList:
        searchList = []
        for x in itemList:
            if x not in searchList:
                searchList.append(x)
            if x == itemList[len(itemList)-1]:
                result.append(searchList)
    return result

#export workflow data
@router.post("/compile")
async def exportData(request: Request, current_user: models.User = Depends(dep.get_current_active_user)):
    try:
        payload_as_json = await request.json()
        inputCon_list = []
        outputCon_list = []
        # print(payload_as_json)
        for key, val in payload_as_json.items():
            for val_key, val_val in val.items():
                if val_key == "inputs":
                    for I_val in val_val.values():
                        inputCon_list.append([key, I_val['connections'][0]['node']])
                elif val_key == "outputs":
                    for O_val in val_val.values():
                        if len(O_val['connections']) != 0:
                            outputCon_list.append([key, O_val['connections'][0]['node']])
                        
        nodeObject_list = linkList(outputCon_list)
        # print(nodeObject_list)
        for item in nodeObject_list:
            item_file = payload_as_json[item[0]]["data"]["file"].replace('C:\\fakepath\\', '').replace('.csv', '')
            lastNode = payload_as_json[item[len(item)-1]]["name"].replace(' ', '')
            print(item_file, lastNode)
            with open(f"workflow/data/{current_user.username}_{item_file}.txt", 'w') as f:
                f.write(item_file)
            process = Popen(['snakemake',f'workflow/data/{lastNode}_{current_user.username}_{item_file}.csv','-j'], stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
        message = "success"
    except json.JSONDecodeError:
        payload_as_json = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": nodeObject_list}

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

    if FILE_NAME.find("ScatterPlot") != -1:
        img = Image.open(FILE_PATH)
        img_converted = from_image_to_bytes(img)
        return JSONResponse(img_converted)
    else:
        return FileResponse(FILE_PATH)

    # 차후 개발 방향
    # workflow DB에서 가장 최근에 생성된 Column 가져옴
    # 노드 정보들을 통해 file_list 안에 해당 노드 결과들이 생성되었는지 파악
