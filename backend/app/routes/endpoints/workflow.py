from fastapi import APIRouter, Body, Depends, HTTPException, Request
from fastapi.encoders import jsonable_encoder

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