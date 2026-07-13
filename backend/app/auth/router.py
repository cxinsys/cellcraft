from typing import Any

from fastapi import APIRouter, Depends
from fastapi.security import OAuth2PasswordRequestForm

from sqlalchemy.orm import Session

from app.auth import schemas as model
from app import models
from app.user import schemas as user
from app.auth import service as auth_service
from app.user import service as user_service
from app.auth import deps as dep

router = APIRouter()


# create New User
@router.post("/register", response_model=user.UserProfile)
def create_user(
    *,
    db: Session = Depends(dep.get_db),
    user_in: user.UserCreate,
) -> Any:
    return user_service.register(db=db, user_in=user_in)


# Login + JWT 발급
@router.post("/login/access-token", response_model=model.Token)
def login_access_token(
    db: Session = Depends(dep.get_db), form_data: OAuth2PasswordRequestForm = Depends()
) -> Any:
    return auth_service.login(db=db, form_data=form_data)


# read Current User
@router.get("/me", response_model=user.UserProfile)
def read_user_me(
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    return user_service.read_me(current_user=current_user)


# update user plugins
@router.post("/plugins", response_model=user.UserProfile)
def update_user_plugins(
    *,
    db: Session = Depends(dep.get_db),
    user_in: user.UserUpdate,
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    return user_service.update_plugins(db=db, current_user=current_user, user_in=user_in)
