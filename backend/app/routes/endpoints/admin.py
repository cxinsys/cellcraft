from fastapi import APIRouter, Depends, HTTPException
from typing import List,Optional,Any
from sqlalchemy.orm import Session
import psutil

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
    # conditions: Conditions, # 이거 안됨
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
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

@router.get("/system/stats", response_model=Any)
def get_system_stats():
    # CPU 사용량
    cpu_usage = psutil.cpu_percent(interval=1)
    
    # 메모리 사용량
    memory_info = psutil.virtual_memory()
    memory_usage = memory_info.percent
    total_memory = memory_info.total
    available_memory = memory_info.available
    used_memory = memory_info.used

    # 반환할 데이터
    system_stats = {
        "cpu_usage_percent": cpu_usage,
        "memory_usage_percent": memory_usage,
        "total_memory_bytes": total_memory,
        "available_memory_bytes": available_memory,
        "used_memory_bytes": used_memory
    }

    return system_stats