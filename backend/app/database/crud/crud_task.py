from datetime import datetime
from sqlalchemy.orm import Session, joinedload
from sqlalchemy.orm.exc import MultipleResultsFound
from app.database import models
from app.database.conn import get_db_session
from app.database.crud import crud_plugin
import logging

logger = logging.getLogger(__name__)

def start_task(user_id: int, task_id: str, workflow_id: int, start_time: datetime, algorithm_id: str = None, plugin_name: str = None, task_type: str = None, plugin_image_uri: str = None, plugin_id: int = None):
    """
    Start a new task with optional plugin_id for enhanced plugin tracking.

    Args:
        user_id: User ID
        task_id: Celery task ID
        workflow_id: Workflow ID
        start_time: Task start time
        algorithm_id: Algorithm ID (optional)
        plugin_name: Plugin name (kept for backward compatibility)
        task_type: Task type (compile, visualization, etc.)
        plugin_image_uri: Plugin image URI
        plugin_id: Plugin ID for foreign key relationship (optional)
    """
    with get_db_session() as db:
        # If plugin_id is not provided but plugin_name is, try to find the plugin_id
        if plugin_id is None and plugin_name:
            plugin_id = find_plugin_id_by_name(db, plugin_name)

        db_task = models.Task(
            user_id=user_id,
            task_id=task_id,
            workflow_id=workflow_id,
            algorithm_id=algorithm_id,
            start_time=start_time,
            status='RUNNING',
            plugin_name=plugin_name,
            task_type=task_type,
            plugin_image_uri=plugin_image_uri,
            plugin_id=plugin_id
        )
        db.add(db_task)
        # commit은 context manager가 자동 처리

        if plugin_id:
            logger.info(f"Task {task_id} started with plugin_id {plugin_id}")
        else:
            logger.info(f"Task {task_id} started without plugin_id (legacy mode)")

def end_task(user_id: int, task_id: str, end_time: datetime, status: str):
    with get_db_session() as db:
        try:
            task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).one()
            task.end_time = end_time
            task.status = status
            # commit은 context manager가 자동 처리
        except MultipleResultsFound:
            # 여러 개의 결과가 발견되었을 때의 처리를 합니다.
            # 예를 들어, 첫 번째 결과를 사용하거나, 오류 메시지를 출력하거나, 데이터를 정리하는 등의 작업을 수행할 수 있습니다.
            task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).first()
            # 첫 번째 결과를 사용하기로 결정했으면 위와 같이 .first()를 사용하여 첫 번째 결과를 가져올 수 있습니다.
            # Note: 원래 코드에서도 MultipleResultsFound 케이스에서는 업데이트하지 않음 (데이터 품질 이슈로 별도 처리 필요)

def get_user_task(db: Session, id: int):
    return db.query(models.Task).filter(models.Task.user_id == id).all()

def get_user_task_with_plugin(db: Session, user_id: int):
    """
    Get user tasks with joined plugin information and workflow data for enhanced API responses.
    Uses eager loading to prevent N+1 query problems.
    
    Args:
        db: Database session
        user_id: User ID
        
    Returns:
        List of Task objects with joined Plugin and Workflow data
    """
    return db.query(models.Task).filter(
        models.Task.user_id == user_id
    ).options(
        joinedload(models.Task.plugin),
        joinedload(models.Task.workflows)
    ).all()

def get_task_by_task_id(db: Session, task_id: str):
    return db.query(models.Task).filter(models.Task.task_id == task_id).first()

def delete_user_task(db: Session, user_id: int, task_id: str):
    target_task = db.query(models.Task).filter(models.Task.task_id == task_id, models.Task.user_id == user_id).first()
    if target_task is None:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail=f"Task {task_id} not found")
    db.delete(target_task)
    db.commit()
    return target_task

def record_plugin_image_uri(task_id: str, plugin_image_uri: str, user_id: int = None):
    """
    Record plugin image URI for an existing task.

    Args:
        task_id: Task ID to update
        plugin_image_uri: The plugin image URI used for execution
        user_id: Optional user ID for additional filtering
    """
    with get_db_session() as db:
        # Build query filters
        filters = [models.Task.task_id == task_id]
        if user_id:
            filters.append(models.Task.user_id == user_id)

        task = db.query(models.Task).filter(*filters).first()
        if task:
            task.plugin_image_uri = plugin_image_uri
            # commit은 context manager가 자동 처리
            logger.info(f"Updated plugin_image_uri for task {task_id}: {plugin_image_uri}")
        else:
            logger.warning(f"Task not found for plugin_image_uri update: {task_id}")

def find_plugin_id_by_name(db: Session, plugin_name: str) -> int:
    """
    Helper function to find plugin_id by plugin name for backward compatibility.
    Prioritizes official plugins over local ones.
    
    Args:
        db: Database session
        plugin_name: Plugin name to search for
        
    Returns:
        Plugin ID if found, None otherwise
    """
    try:
        # Try to find official plugin first (prioritize official over local)
        plugin = db.query(models.Plugin).filter(
            models.Plugin.name == plugin_name,
            models.Plugin.source == "official"
        ).first()
        
        # If no official plugin found, try local plugins
        if not plugin:
            plugin = db.query(models.Plugin).filter(
                models.Plugin.name == plugin_name,
                models.Plugin.source == "local"
            ).first()
            
        # If still no plugin found, try any plugin with that name
        if not plugin:
            plugin = db.query(models.Plugin).filter(
                models.Plugin.name == plugin_name
            ).first()
            
        return plugin.id if plugin else None
        
    except Exception as e:
        logger.warning(f"Failed to find plugin_id for plugin_name {plugin_name}: {e}")
        return None
