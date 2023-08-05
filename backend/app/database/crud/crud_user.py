from typing import Any, Dict, Optional, Union, List

from app.database.schemas.admin import Conditions  # Conditions를 정의한 python 파일을 import

from sqlalchemy.orm import Session
from sqlalchemy import asc, desc, or_
from app.database import models
from app.database.schemas import user
from app.common.security import get_password_hash, verify_password


def get_user(db: Session, id: int):
    return db.query(models.User).filter(models.User.id == id).first()

def get_filtered_users(db: Session, conditions: Conditions) -> List[models.User]:
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    sort_column = getattr(models.User, sort, None)
    if not sort_column:
        raise ValueError(f"Sort column {sort} does not exist on User model")

    query = db.query(models.User)
    
    if searchTerm:
        query = query.filter(
            or_(
                models.User.username.like(f"%{searchTerm}%"),
                models.User.email.like(f"%{searchTerm}%")
            )
        )
    
    order_func = asc if order == 'asc' else desc
    return query.order_by(order_func(sort_column)).offset(amount * (page_num - 1)).limit(amount).all()

def get_users_count(db: Session) -> int:
    return db.query(models.User).count()

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