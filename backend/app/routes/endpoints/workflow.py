from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File
from fastapi.encoders import jsonable_encoder
from typing import List, Union
import os

from json import JSONDecodeError

router = APIRouter()

#export workflow data
@router.post("/export")
async def exportData(request: Request):
    try:
        payload_as_json = await request.json()
        message = "success"
    except JSONDecodeError:
        payload_as_json = None
        message = "Received data is not a valid JSON"
    return {"message": message, "recived_data": payload_as_json}

#workflow file-upload
@router.post("/upload")
async def fileUpload(files: List[UploadFile] = File()):
    UPLOAD_DIRECTORY = './'
    for item_file in files:
        contents = await item_file.read()
        with open(os.path.join(UPLOAD_DIRECTORY, item_file.filename), "wb") as f:
            f.write(contents)
    return {"filename" : [item_file.filename for item_file in files]}
