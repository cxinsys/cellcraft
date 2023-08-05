from fastapi import APIRouter, Depends, HTTPException
from typing import List,Optional,Any
from sqlalchemy.orm import Session

from app.routes import dep  # get_db를 정의한 python 파일을 import
from app.database.schemas.user import User  # User를 정의한 python 파일을 import
from app.database.schemas.admin import Conditions  # Conditions를 정의한 python 파일을 import
from app.database.crud import crud_user  # get_filtered_users를 정의한 python 파일을 import

from app.database import models

router = APIRouter()

@router.get("/users", response_model=List[User])
def get_filtered_users(
    *,
    db: Session = Depends(dep.get_db),
    amount: int = 20,
    page_num: int = 1,
    sort: str = "id",
    order: str = "asc",
    searchTerm: Optional[str] = None
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    users = crud_user.get_filtered_users(db, conditions)

    if users is None:
        raise HTTPException(status_code=404, detail="Users not found")
    return users

@router.get("/users_count", response_model=Any)
def get_users_count(
    *,
    db: Session = Depends(dep.get_db),
):
    users_num = crud_user.get_users_count(db)
    return users_num