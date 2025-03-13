from sqlalchemy.orm import Session
from sqlalchemy import asc, desc, or_
from typing import List
from app.database import models

def get_filtered_users(db: Session, conditions) -> List[models.User]:
    return get_filtered_data(db, models.User, conditions, ['username', 'email'])

def get_filtered_files(db: Session, conditions) -> List[models.File]:
    return get_filtered_data(db, models.File, conditions, ['file_name', 'folder'])

def get_filtered_workflows(db: Session, conditions) -> List[models.Workflow]:
    return get_filtered_data(db, models.Workflow, conditions, ['title'])

def get_filtered_tasks(db: Session, conditions) -> List[models.Task]:
    return get_filtered_data(db, models.Task, conditions, ['task_id', 'status'])

def get_filtered_plugins(db: Session, conditions) -> List[models.Plugin]:
    return get_filtered_data(db, models.Plugin, conditions, ['name', 'author'])

def get_filtered_data(db: Session, model, conditions, searchable_fields):
    """
    일반화된 데이터 필터링 함수
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    sort_column = getattr(model, sort, None)
    if not sort_column:
        raise ValueError(f"Sort column {sort} does not exist on {model.__name__} model")

    query = db.query(model)

    if searchTerm:
        search_filters = [getattr(model, field).like(f"%{searchTerm}%") for field in searchable_fields]
        query = query.filter(or_(*search_filters))

    order_func = asc if order == 'asc' else desc
    return query.order_by(order_func(sort_column)).offset(amount * (page_num - 1)).limit(amount).all()

def get_users_count(db: Session) -> int:
    return get_count(db, models.User)

def get_files_count(db: Session) -> int:
    return get_count(db, models.File)

def get_workflows_count(db: Session) -> int:
    return get_count(db, models.Workflow)

def get_tasks_count(db: Session) -> int:
    return get_count(db, models.Task)

def get_plugins_count(db: Session) -> int:
    return get_count(db, models.Plugin)

def get_count(db: Session, model) -> int:
    """
    주어진 모델의 총 개수를 반환하는 공통 함수
    """
    return db.query(model).count()