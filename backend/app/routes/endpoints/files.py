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
from app.database.schemas.file import FileCreate, FileDelete, FileUpdate, FileFind, FolderFind


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
    UPLOAD_DIRECTORY = f'./user/{current_user.username}/data'
    for item_file in files:
        contents = await item_file.read()
        print(item_file.filename.split("_"))
        folder_file = item_file.filename.split('_')
        with open(os.path.join(UPLOAD_DIRECTORY, folder_file[1]), "wb") as f:
            f.write(contents)
        #file 있는지 여부 검증
        #file 있으면 에러
        #file create
        #upload to DB
        user_file = crud_file.get_user_file(db, current_user.id, folder_file[1])
        if user_file:
            raise HTTPException(
                status_code=400,
                detail="this file already exists in your files",
                )
        print(len(contents))
        create_file = await crud_file.create_file(db, folder_file[1], len(contents), UPLOAD_DIRECTORY, folder_file[0], current_user.id)
        
    return create_file

#User Files get
@router.get("/me")
def read_user_folder(
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any:
    USER_FOLDER = f'./user/{current_user.username}/data'
    res = []
    for (dir_path, dir_names, file_names) in os.walk(USER_FOLDER):
        folder = dir_path.replace(USER_FOLDER , 'data')
        res.extend({folder : file_names}.items())
        print(f"Directories: {dir_path}, Files: {file_names}")
    print(res)

    # print(user_files)
    return res

#User File find
@router.post("/find")
def find_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        return user_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )

#User find File of folder
@router.post("/folder")
def find_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    folder: FolderFind,
    ) -> Any:
    user_file = crud_file.get_user_folder(db, current_user.id, folder.folder_name)
    if user_file:
        return user_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this folder not exists in your folders",
                )

#User File delete
@router.post("/delete")
def delete_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileDelete,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        delete_file = crud_file.delete_user_file(db, current_user.id, fileInfo.file_name)
        print(delete_file)
        return delete_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )

#User File Update
@router.post("/update")
def update_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileUpdate,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        update_file = crud_file.update_user_file(db, current_user.id, fileInfo)
        print(update_file)
        return update_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )