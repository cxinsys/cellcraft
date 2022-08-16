from typing import Any, Dict, Optional, Union

from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import user
from app.common.security import get_password_hash, verify_password


def get_user(db: Session, id: int):
    return db.query(models.User).filter(models.User.id == id).first()


def get_user_by_email(db: Session, email: str):
    return db.query(models.User).filter(models.User.email == email).first()


def create_user(db: Session, user: user.UserCreate) -> models.User:
    db_user = models.User(
        email=user.email, 
        hashed_password=get_password_hash(user.password),
        username=user.username,
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