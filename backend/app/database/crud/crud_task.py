from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import MultipleResultsFound
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
    except MultipleResultsFound:
        # 여러 개의 결과가 발견되었을 때의 처리를 합니다.
        # 예를 들어, 첫 번째 결과를 사용하거나, 오류 메시지를 출력하거나, 데이터를 정리하는 등의 작업을 수행할 수 있습니다.
        task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).first()
        # 첫 번째 결과를 사용하기로 결정했으면 위와 같이 .first()를 사용하여 첫 번째 결과를 가져올 수 있습니다.
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


