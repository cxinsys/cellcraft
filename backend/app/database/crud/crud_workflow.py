from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import workflow

def create_workflow(db: Session, title: str, workflow_info: dict, nodes: list, linked_nodes: list, user_id: int) -> models.Workflow:
    db_workflow = models.Workflow(
        title = title,
        workflow_info = workflow_info,
        nodes = nodes,
        linked_nodes = linked_nodes,
        user_id = user_id
        )
    db.add(db_workflow)
    db.commit()
    db.refresh(db_workflow)
    return db_workflow

def get_user_workflows(db: Session, user_id: int):
    return db.query(models.Workflow).filter(models.Workflow.user_id == user_id).all()

def get_user_workflow(db: Session, user_id: int, id: int):
    return db.query(models.Workflow).filter(models.Workflow.id == id, models.Workflow.user_id == user_id).first()