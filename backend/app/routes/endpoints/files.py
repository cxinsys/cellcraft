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

router = APIRouter()

#workflow file-upload
@router.post("/upload")
def fileUpload(
    *,
    db: Session = Depends(dep.get_db),
    files: List[UploadFile] = File(),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    #upload to directory
    UPLOAD_DIRECTORY = './'
    for item_file in files:
        contents = item_file.read()
        with open(os.path.join(UPLOAD_DIRECTORY, item_file.filename), "wb") as f:
            f.write(contents)
        
    #     file_name = Column(String, nullable=False)
    # file_size = Column(String, nullable=False)
    # file_path = Column(String, nullable=False)
    # user_id = Column(Integer, ForeignKey("users.id"))
        print(len(contents))
        #file 있는지 여부 검증
        #file 있으면 에러
        #file create
        
        # user_file = crud_file.get_user_file(db, current_user.user_id, item_file.filename)
    #upload to DB
    
    return {"filename" : [item_file.filename for item_file in files]}