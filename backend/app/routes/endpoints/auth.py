from typing import Any, List
from datetime import timedelta

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from fastapi.security import OAuth2PasswordRequestForm

from pydantic.networks import EmailStr
from sqlalchemy.orm import Session

import os

from app import model
from app.database import models
from app.database.crud import crud_user
from app.database.schemas import user
from app.common.security import create_access_token
from app.common.config import settings
from app.routes import dep

router = APIRouter()

#create New User
@router.post("/register", response_model=user.User)
def create_user(
    *,
    db: Session = Depends(dep.get_db),
    user_in: user.UserCreate,
) -> Any:
    try:
        user = crud_user.get_user_by_email(db, email=user_in.email)
        if user:
            raise HTTPException(
                status_code=400,
                detail="this email already exists in the system",
            )
        user = crud_user.create_user(db, user=user_in)
        USER_DIRECTORY_NAME = './user/' + user_in.username + '/data'
        os.makedirs(USER_DIRECTORY_NAME)
        #회원가입 시 보내는 확인 이메일
        # if user_in.email:
        #     send_new_account_email(
        #         email_to=user_in.email, username=user_in.email, password=user_in.password
        #     )
        return user
    except FileExistsError as err:
        return err

#Login + JWT 발급
@router.post("/login/access-token", response_model=model.Token)
def login_access_token(
    db: Session = Depends(dep.get_db), form_data: OAuth2PasswordRequestForm = Depends()
) -> Any:
    user = crud_user.authenticate(
        db, email=form_data.username, password=form_data.password
    )
    if not user:
        raise HTTPException(status_code=400, detail="Incorrect email or password")
    elif not crud_user.is_active(user):
        raise HTTPException(status_code=400, detail="please login to activate")
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    return {
        "access_token": create_access_token(
            user.id, expires_delta=access_token_expires
        ),
        "token_type": "bearer",
        "user_info": {
            "is_superuser": crud_user.is_superuser(user),
            "email": user.email,
            "is_active": user.is_active
        }
    }

#read Current User
@router.get("/me", response_model=user.UserProfile)
def read_user_me(
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:

    return { "email" : current_user.email, "username" : current_user.username, "is_superuser" : current_user.is_superuser}

#update user plugins
@router.post("/plugins", response_model=user.User)
def update_user_plugins(
    *,
    db: Session = Depends(dep.get_db),
    user_in: user.UserUpdate,
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    user = crud_user.update_user(db, user_id=current_user.id, user=user_in)
    return user