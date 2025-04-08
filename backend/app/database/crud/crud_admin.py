from sqlalchemy.orm import Session
from sqlalchemy import asc, desc, or_, text
from typing import List, Tuple
from app.database import models
from app.database.schemas.admin import Conditions

def get_filtered_users(db: Session, conditions: Conditions) -> Tuple[List[models.User], int]:
    """
    User 목록을 필터링, 정렬, 페이지네이션하여 반환
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    # 정렬 컬럼 매핑
    sort_mapping = {
        'username': text('users.username'),
        'email': text('users.email'),
        'id': text('users.id')
    }
    
    # 기본 정렬 컬럼 설정
    sort_column = sort_mapping.get(sort, text('users.id'))
    sort_order = asc if order == 'asc' else desc

    # 기본 쿼리 생성
    query = db.query(models.User)

    # 검색 조건 적용
    if searchTerm:
        query = query.filter(
            or_(
                models.User.username.ilike(f'%{searchTerm}%'),
                models.User.email.ilike(f'%{searchTerm}%')
            )
        )

    # 전체 개수 계산
    total_count = query.count()

    # 정렬 및 페이지네이션 적용
    users = query.order_by(sort_order(sort_column)) \
                .offset((page_num - 1) * amount) \
                .limit(amount) \
                .all()

    return users, total_count

def get_filtered_files(db: Session, conditions: Conditions) -> Tuple[List[models.File], int]:
    """
    File 목록을 필터링, 정렬, 페이지네이션하여 반환
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    # 정렬 컬럼 매핑
    sort_mapping = {
        'fileName': text('files.file_name'),
        'folder': text('files.folder'),
        'id': text('files.id'),
        'username': text('users.username'),
        'fileSize': text('files.file_size')
    }
    
    # 기본 정렬 컬럼 설정
    sort_column = sort_mapping.get(sort, text('files.id'))
    sort_order = asc if order == 'asc' else desc

    # 기본 쿼리 생성 (User와 join)
    query = db.query(
        models.File,
        models.User.username.label('username')
    ).join(
        models.User, models.File.user_id == models.User.id
    )

    # 검색 조건 적용
    if searchTerm:
        query = query.filter(
            or_(
                models.File.file_name.ilike(f'%{searchTerm}%'),
                models.File.folder.ilike(f'%{searchTerm}%'),
                models.User.username.ilike(f'%{searchTerm}%')
            )
        )

    # 전체 개수 계산
    total_count = query.count()

    # 정렬 및 페이지네이션 적용
    files = query.order_by(sort_order(sort_column)) \
                .offset((page_num - 1) * amount) \
                .limit(amount) \
                .all()

    return files, total_count

def get_filtered_workflows(db: Session, conditions: Conditions) -> Tuple[List[models.Workflow], int]:
    """
    Workflow 목록을 필터링, 정렬, 페이지네이션하여 반환
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    # 정렬 컬럼 매핑
    sort_mapping = {
        'title': text('workflows.title'),
        'id': text('workflows.id'),
        'username': text('users.username'),
        'updated_at': text('workflows.updated_at')
    }
    
    # 기본 정렬 컬럼 설정
    sort_column = sort_mapping.get(sort, text('workflows.id'))
    sort_order = asc if order == 'asc' else desc

    # 기본 쿼리 생성 (User와 join)
    query = db.query(
        models.Workflow,
        models.User.username.label('username')
    ).join(
        models.User, models.Workflow.user_id == models.User.id
    )

    # 검색 조건 적용
    if searchTerm:
        query = query.filter(
            or_(
                models.Workflow.title.ilike(f'%{searchTerm}%'),
                models.User.username.ilike(f'%{searchTerm}%')
            )
        )

    # 전체 개수 계산
    total_count = query.count()

    # 정렬 및 페이지네이션 적용
    workflows = query.order_by(sort_order(sort_column)) \
                    .offset((page_num - 1) * amount) \
                    .limit(amount) \
                    .all()

    return workflows, total_count

def get_filtered_plugins(db: Session, conditions: Conditions) -> Tuple[List[models.Plugin], int]:
    """
    Plugin 목록을 필터링, 정렬, 페이지네이션하여 반환
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    # 정렬 컬럼 매핑
    sort_mapping = {
        'name': text('plugins.name'),
        'author': text('plugins.author'),
        'id': text('plugins.id')
    }
    
    # 기본 정렬 컬럼 설정
    sort_column = sort_mapping.get(sort, text('plugins.id'))
    sort_order = asc if order == 'asc' else desc

    # 기본 쿼리 생성
    query = db.query(models.Plugin)

    # 검색 조건 적용
    if searchTerm:
        query = query.filter(
            or_(
                models.Plugin.name.ilike(f'%{searchTerm}%'),
                models.Plugin.author.ilike(f'%{searchTerm}%')
            )
        )

    # 전체 개수 계산
    total_count = query.count()

    # 정렬 및 페이지네이션 적용
    plugins = query.order_by(sort_order(sort_column)) \
                  .offset((page_num - 1) * amount) \
                  .limit(amount) \
                  .all()

    return plugins, total_count

def get_users_count(db: Session) -> int:
    return db.query(models.User).count()

def get_files_count(db: Session) -> int:
    return db.query(models.File).count()

def get_workflows_count(db: Session) -> int:
    return db.query(models.Workflow).count()

def get_plugins_count(db: Session) -> int:
    return db.query(models.Plugin).count()

def get_filtered_tasks(db: Session, conditions: Conditions) -> Tuple[List[models.Task], int]:
    """
    Task 목록을 필터링, 정렬, 페이지네이션하여 반환
    """
    amount = conditions.amount
    page_num = conditions.page_num
    sort = conditions.sort
    order = conditions.order
    searchTerm = conditions.searchTerm

    # 정렬 컬럼 매핑
    sort_mapping = {
        'username': text('users.username'),
        'workflowTitle': text('workflows.title'),
        'status': text('tasks.status'),
        'time': text('tasks.start_time'),
        'id': text('tasks.id')
    }
    
    # 기본 정렬 컬럼 설정
    sort_column = sort_mapping.get(sort, text('tasks.id'))
    sort_order = asc if order == 'asc' else desc

    # 기본 쿼리 생성
    query = db.query(
        models.Task,
        models.User.username.label('username'),
        models.Workflow.title.label('workflow_title')
    ).join(
        models.User, models.Task.user_id == models.User.id
    ).join(
        models.Workflow, models.Task.workflow_id == models.Workflow.id
    )

    # 검색 조건 적용
    if searchTerm:
        query = query.filter(
            or_(
                models.User.username.ilike(f'%{searchTerm}%'),
                models.Workflow.title.ilike(f'%{searchTerm}%')
            )
        )

    # 전체 개수 계산
    total_count = query.count()

    # 정렬 및 페이지네이션 적용
    tasks = query.order_by(sort_order(sort_column)) \
                .offset((page_num - 1) * amount) \
                .limit(amount) \
                .all()

    return tasks, total_count

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

def get_tasks_count(db: Session) -> int:
    """
    전체 Task 개수를 반환
    """
    return db.query(models.Task).count()

def update_user(db: Session, user_id: int, user_data: dict) -> models.User:
    """
    사용자 정보를 업데이트
    """
    user = db.query(models.User).filter(models.User.id == user_id).first()
    if not user:
        return None
    
    for key, value in user_data.items():
        if hasattr(user, key):
            setattr(user, key, value)
    
    db.commit()
    db.refresh(user)
    return user

def delete_user(db: Session, user_id: int) -> bool:
    """
    사용자를 삭제
    """
    user = db.query(models.User).filter(models.User.id == user_id).first()
    if not user:
        return False
    
    db.delete(user)
    db.commit()
    return True

def delete_file(db: Session, file_id: int) -> bool:
    """
    파일을 삭제
    """
    file = db.query(models.File).filter(models.File.id == file_id).first()
    if not file:
        return False
    
    db.delete(file)
    db.commit()
    return True

def delete_workflow(db: Session, workflow_id: int) -> bool:
    """
    워크플로우를 삭제
    """
    workflow = db.query(models.Workflow).filter(models.Workflow.id == workflow_id).first()
    if not workflow:
        return False
    
    db.delete(workflow)
    db.commit()
    return True

def cancel_task(db: Session, task_id: int) -> bool:
    """
    태스크를 취소
    """
    task = db.query(models.Task).filter(models.Task.id == task_id).first()
    if not task:
        return False
    
    task.status = "cancelled"
    db.commit()
    return True

def install_plugin_dependencies(db: Session, plugin_id: int) -> bool:
    """
    플러그인 의존성을 설치
    """
    plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
    if not plugin:
        return False
    
    # TODO: 실제 의존성 설치 로직 구현
    # 예: subprocess를 사용하여 pip install 실행
    return True