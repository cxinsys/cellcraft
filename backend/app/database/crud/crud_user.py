from typing import Any, Dict, Optional, Union, List

from app.database.schemas.admin import Conditions  # Conditions를 정의한 python 파일을 import

from sqlalchemy.orm import Session
from sqlalchemy import asc, desc, or_
from app.database import models
from app.database.schemas import user
from app.common.security import get_password_hash, verify_password

def get_user(db: Session, id: int):
    return db.query(models.User).filter(models.User.id == id).first()

def get_users_count(db: Session) -> int:
    return db.query(models.User).count()

def get_user_by_email(db: Session, email: str):
    return db.query(models.User).filter(models.User.email == email).first()

def update_user(db: Session, user_id: int, user: user.UserUpdate) -> models.User:
    db_user = db.query(models.User).filter(models.User.id == user_id).first()
    if user.password:
        db_user.hashed_password = get_password_hash(user.password)
    db.commit()
    db.refresh(db_user)
    return db_user

def create_user(db: Session, user: user.UserCreate) -> models.User:
    db_user = models.User(
        email=user.email, 
        hashed_password=get_password_hash(user.password),
        username=user.username,
        )
    db.add(db_user)
    db.commit()
    db.refresh(db_user)

    # 모든 플러그인을 조회하여 사용자와 연결
    all_plugins = db.query(models.Plugin).all()
    for plugin in all_plugins:
        db_user.plugins.append(plugin)
    
    db.commit()
    db.refresh(db_user)
    return db_user

def create_superuser(db: Session, user: user.UserCreate) -> models.User:
    db_user = models.User(
        email=user.email, 
        hashed_password=get_password_hash(user.password),
        username=user.username,
        is_superuser=True,
        is_active=True  # 여기에 is_active 값을 추가합니다.
        )
    db.add(db_user)
    db.commit()
    db.refresh(db_user)
    return db_user

def authenticate(db: Session, *, email: str, password: str) -> Optional[models.User]:
        user = get_user_by_email(db, email=email)
        if not user:
            print('not email')
            return None
        if not verify_password(password, user.hashed_password):
            print('not password')
            return None
        return user

def is_active(user: models.User) -> bool:
        return user.is_active

def is_superuser(user: models.User) -> bool:
        return user.is_superuser