from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import workflow

def create_workflow(db: Session, title: str, thumbnail:str, workflow_info: dict, user_id: int) -> models.Workflow:
    db_workflow = models.Workflow(
        title = title,
        thumbnail = thumbnail,
        workflow_info = workflow_info,
        user_id = user_id
        )
    db.add(db_workflow)
    db.commit()
    db.refresh(db_workflow)
    return db_workflow

def update_workflow(db: Session, user_id: int, workflow_id: int, title: str = None, thumbnail: str = None, workflow_info: dict = None) -> models.Workflow:
    target_workflow = db.query(models.Workflow).filter(models.Workflow.id == workflow_id, models.Workflow.user_id == user_id).first()
    if title:
        target_workflow.title = title
    if thumbnail:
        target_workflow.thumbnail = thumbnail
    if workflow_info:
        target_workflow.workflow_info = workflow_info
    db.commit()
    db.refresh(target_workflow)
    return target_workflow

def get_user_workflows(db: Session, user_id: int):
    return db.query(models.Workflow).filter(models.Workflow.user_id == user_id).all()

def get_user_workflow(db: Session, user_id: int, id: int):
    return db.query(models.Workflow).filter(models.Workflow.id == id, models.Workflow.user_id == user_id).first()

def delete_user_workflow(db: Session, user_id: int, workflow_id: id):
    target_workflow = db.query(models.Workflow).filter(models.Workflow.id == workflow_id, models.Workflow.user_id == user_id).first()
    db.delete(target_workflow)
    db.commit()
    return target_workflow

def get_workflow_by_id(db: Session, workflow_id: int) -> models.Workflow:
    return db.query(models.Workflow).filter(models.Workflow.id == workflow_id).first()