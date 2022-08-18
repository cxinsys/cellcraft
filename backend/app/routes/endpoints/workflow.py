from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File, responses
from fastapi.responses import FileResponse
from fastapi.encoders import jsonable_encoder
from typing import List, Union
from subprocess import Popen,PIPE
import os
import json

from app.database.schemas.workflow import WorkflowResult

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
async def exportData(request: Request):
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
            with open(f"workflow/data/build_{item_file}.txt", 'w') as f:
                f.write(item_file)
            process = Popen(['snakemake',f'workflow/data/{lastNode}_{item_file}.csv','-j'], stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
        message = "success"
    except json.JSONDecodeError:
        payload_as_json = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": nodeObject_list}

@router.post("/result")
def checkResult(filename: WorkflowResult):
    PATH_COMPILE_RESULT = './workflow/result'
    file_list = os.listdir(PATH_COMPILE_RESULT)
    print(file_list)
    FILE_NAME = filename.filename
    print(FILE_NAME)
    for item_file in file_list:
        if FILE_NAME in item_file:
            FILE_NAME = item_file
    
    FILE_PATH = PATH_COMPILE_RESULT + '/' + FILE_NAME
    
    # workflow DB에서 가장 최근에 생성된 Column 가져옴
    # 노드 정보들을 통해 file_list 안에 해당 노드 결과들이 생성되었는지 파악
    return FileResponse(FILE_PATH)
