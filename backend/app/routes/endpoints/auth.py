from typing import Any, List
from datetime import timedelta

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from fastapi.security import OAuth2PasswordRequestForm

from pydantic.networks import EmailStr
from sqlalchemy.orm import Session

from app import model
from app.database import crud, models, schemas
from app.common.security import create_access_token
from app.common.config import settings
from app.routes import dep

router = APIRouter()

#create New User
@router.post("/register", response_model=schemas.User)
def create_user(
    *,
    db: Session = Depends(dep.get_db),
    user_in: schemas.UserCreate,
) -> Any:
    user = crud.get_user_by_email(db, email=user_in.email)
    if user:
        raise HTTPException(
            status_code=400,
            detail="this email already exists in the system",
        )
    user = crud.create_user(db, user=user_in)
    #회원가입 시 보내는 확인 이메일
    # if user_in.email:
    #     send_new_account_email(
    #         email_to=user_in.email, username=user_in.email, password=user_in.password
    #     )
    return user

#Login + JWT 발급
@router.post("/login/access-token", response_model=model.Token)
def login_access_token(
    db: Session = Depends(dep.get_db), form_data: OAuth2PasswordRequestForm = Depends()
) -> Any:
    user = crud.authenticate(
        db, email=form_data.username, password=form_data.password
    )
    if not user:
        raise HTTPException(status_code=400, detail="Incorrect email or password")
    elif not crud.is_active(user):
        raise HTTPException(status_code=400, detail="please login to activate")
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    return {
        "access_token": create_access_token(
            user.id, expires_delta=access_token_expires
        ),
        "token_type": "bearer",
    }

#read Current User
@router.get("/me", response_model=schemas.User)
def read_user_me(
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:

    return current_user