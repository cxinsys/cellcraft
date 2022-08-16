from typing import Any
from venv import create
from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File
from fastapi.encoders import jsonable_encoder
from typing import List, Union
from sqlalchemy.orm import Session
import os
from json import JSONDecodeError

from app.routes import dep
from app.database.crud import crud_file
from app.database import models
from app.database.schemas.file import FileCreate

router = APIRouter()

#workflow file-upload
@router.post("/upload")
async def fileUpload(
    *,
    db: Session = Depends(dep.get_db),
    files: List[UploadFile] = File(),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    #upload to directory
    UPLOAD_DIRECTORY = './'
    for item_file in files:
        contents = await item_file.read()
        with open(os.path.join(UPLOAD_DIRECTORY, item_file.filename), "wb") as f:
            f.write(contents)
        #file 있는지 여부 검증
        #file 있으면 에러
        #file create
        #upload to DB
        user_file = crud_file.get_user_file(db, current_user.id, item_file.filename)
        if user_file:
            raise HTTPException(
                status_code=400,
                detail="this file already exists in your files",
                )
        print(len(contents))
        create_file = await crud_file.create_file(db, item_file, len(contents), current_user.id)
        
    return create_file

#User Files get
@router.get("/me")
def read_user_files(
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any:
    user_files = crud_file.get_user_files(db, current_user.id)
    print(user_files)
    return user_files
