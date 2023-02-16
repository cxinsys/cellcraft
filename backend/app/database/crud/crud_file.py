from typing import Any, Dict, Optional, Union
from fastapi import UploadFile, File
from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import user, file
from app.database.crud import crud_file

def get_user_file(db: Session, id: int, file_name: str):
    return db.query(models.File).filter(models.File.file_name == file_name, models.File.user_id == id).first()

def get_user_folder(db: Session, id: int, folder_name: str):
    return db.query(models.File).filter(models.File.folder == folder_name, models.File.user_id == id).all()

def get_user_files(db: Session, id: int):
    return db.query(models.File).filter(models.File.user_id == id).all()

async def create_file(db: Session, filename: str, filesize: str, filepath: str, folder: str, user_id: int) -> models.File:
    db_file = models.File(
        file_name = filename,
        file_size = filesize,
        file_path = filepath,
        folder = folder,
        user_id = user_id
        )
    db.add(db_file)
    db.commit()
    db.refresh(db_file)
    return db_file

def delete_user_file(db: Session, id: int, file_name: str):
    target_file = db.query(models.File).filter(models.File.file_name == file_name, models.File.user_id == id).first()
    db.delete(target_file)
    db.commit()
    return target_file

def update_user_file(db: Session, id: int, update_info: file.FileUpdate):
    target_file = db.query(models.File).filter(models.File.file_name == update_info.file_name, models.File.user_id == id).first()
    target_file.file_name = update_info.update_name
    db.commit()
    db.refresh(target_file)
    return target_file