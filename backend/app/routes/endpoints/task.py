from fastapi import APIRouter, Depends, HTTPException
from sse_starlette.sse import EventSourceResponse
from typing import Any
from sqlalchemy.orm import Session
import asyncio
from datetime import datetime

from app.common.utils.celery_utils import get_task_info
from app.database.crud import crud_task
from app.routes import dep
from app.database import models

router = APIRouter()                  

@router.get("/info/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    async def event_generator():
        while True:
            if task_id:
                task = get_task_info(task_id)
                if task['task_status'] == 'SUCCESS' or task['task_status'] == 'FAILURE' or task['task_status'] == 'REVOKED' or task['task_status'] == 'RETRY':
                    yield f"{task['task_status']}"
                    break
                print(task['task_status'])
                yield f"{task['task_status']}"
                await asyncio.sleep(5)
            else:
                break
    return EventSourceResponse(event_generator())

@router.get("/monitoring")
async def get_task_monitoring(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any :
    """
    Return the status of the all User Task
    """
    user_task = crud_task.get_user_task(db, current_user.id)
    if user_task:
        return user_task
    else:
        raise HTTPException(
                status_code=400,
                detail="this user not exists task",
                )

@router.delete("/revoke/{task_id}")
def revoke_task(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
    ) -> dict:
    """
    Revoke the task
    """
    from app.main import get_celery_app
    celery = get_celery_app()
    celery.control.revoke(task_id, terminate=True, signal='SIGTERM')
    task = get_task_info(task_id)

    print(task)

    # task_info가 사전 형식인 경우 상태를 접근하는 방식
    task_status = task.get("status")
    if task_status == 'REVOKED':
        return {"message": "Task Revoked", "task_id": task_id}
    else:
        # 태스크 상태를 'REVOKED'로 업데이트
        crud_task.end_task(current_user.id, task_id, datetime.now(), 'REVOKED')
        # 태스크가 업데이트 되었는지 확인
        task = get_task_info(task_id)
        task_status = task.get("status")
        if task_status == 'REVOKED':
            return {"message": "Task Revoked", "task_id": task_id}
        else:
            return {"message": "Task Revoked Failed", "task_id": task_id}

@router.delete("/delete/{task_id}")
def delete_task(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
    ) -> dict:
    """
    Delete the task
    """
    task = crud_task.delete_user_task(db, current_user.id, task_id)
    return {"message": "Task Deleted", "task_id": task_id}
