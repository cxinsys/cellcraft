from datetime import datetime
from sqlalchemy.orm import Session
from app.database import models
from app.database.conn import get_new_engine_and_session

def start_task(user_id: int, task_id: str, workflow_id: int, start_time: datetime):
    db = get_new_engine_and_session()
    try:
        db_task = models.Task(
            user_id=user_id,
            task_id=task_id,
            workflow_id=workflow_id,
            start_time=start_time,
            status='RUNNING'
        )
        db.add(db_task)
        db.commit()
        db.refresh(db_task)
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()

def end_task(user_id: int, task_id: str, end_time: datetime, status: str):
    db = get_new_engine_and_session()
    try:
        task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).one()
        task.end_time = end_time
        task.status = status
        db.commit()
        db.refresh(task)
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()

def get_user_task(db: Session, id: int):
    return db.query(models.Task).filter(models.Task.user_id == id).all()

def delete_user_task(db: Session, user_id: int, task_id: str):
    target_task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).first()
    db.delete(target_task)
    db.commit()
    return target_task


