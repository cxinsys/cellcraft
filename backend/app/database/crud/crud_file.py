from typing import Any, Dict, Optional, Union

from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import user, file
from app.database.crud import crud_file

def get_user_file(db: Session, id: int, file_name: int):
    user_id = db.query(models.User.id).filter(models.User.id == id).first()
    return db.query(models.File).filter(models.File.file_name == file_name, models.File.user_id == user_id)

def create_file(db: Session, file: file.FileCreate, user_id: int) -> models.File:
    db_file = models.File(
        file_name = file.file_name,
        file_size = file.file_size,
        file_path = file.file_path,
        user_id = user_id
        )
    db.add(db_file)
    db.commit()
    db.refresh(db_file)
    return db_file