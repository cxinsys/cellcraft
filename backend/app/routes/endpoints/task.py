from fastapi import APIRouter, Depends, HTTPException
from sse_starlette.sse import EventSourceResponse
from typing import Any, List, Dict
from sqlalchemy.orm import Session
import asyncio
from datetime import datetime

from app.common.utils.celery_utils import get_task_info
from app.common.utils.docker_utils import container_manager
from app.database.crud import crud_task, crud_workflow
from app.routes import dep
from app.database import models
import os
from pathlib import Path

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
    ) -> List[Dict[str, Any]]:
    """
    Return the status of the all User Task with workflow information
    플러그인 빌드 태스크는 제외하고 워크플로우 관련 태스크만 반환
    """
    user_tasks = crud_task.get_user_task(db, current_user.id)
    if not user_tasks:
        raise HTTPException(
            status_code=400,
            detail="this user not exists task",
        )
    
    # 플러그인 빌드 태스크는 제외하고 워크플로우 관련 태스크만 필터링
    workflow_tasks = [task for task in user_tasks if task.task_type != "plugin_build"]
    
    # 각 task에 workflow 정보 및 task_title 추가
    tasks_with_workflow = []
    for task in workflow_tasks:
        workflow = crud_workflow.get_workflow_by_id(db, task.workflow_id)
        
        # task_title 생성 로직
        task_title = None
        if task.plugin_name:
            if task.task_type == 'compile':
                task_title = task.plugin_name
            elif task.task_type == 'visualization':
                task_title = f"{task.plugin_name}-visualization"
            else:
                task_title = task.plugin_name  # 기본값
        
        task_dict = {
            "id": task.id,
            "task_id": task.task_id,  # 실제 Celery task ID 추가
            "workflow_id": task.workflow_id,
            "user_id": task.user_id,
            "status": task.status,
            "start_time": task.start_time,
            "end_time": task.end_time,
            "workflow_title": workflow.title if workflow else None,
            "task_title": task_title
        }
        tasks_with_workflow.append(task_dict)
    
    return tasks_with_workflow

@router.delete("/revoke/{task_id}")
def revoke_task(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
    ) -> dict:
    """
    Revoke the task and cleanup associated containers
    """
    try:
        from app.main import get_celery_app
        celery = get_celery_app()
        
        print(f"Attempting to revoke task {task_id}")
        
        # 1. 먼저 관련 컨테이너 정보 확인
        container_id = container_manager.get_container_id(task_id)
        if container_id:
            print(f"Found associated container {container_id} for task {task_id}")
        
        # 2. Celery 작업 취소
        celery.control.revoke(task_id, terminate=True, signal='SIGTERM')
        print(f"Celery revoke command sent for task {task_id}")
        
        # 3. 컨테이너 매니저를 통한 컨테이너 강제 정리
        container_cleanup_success = False
        if container_id:
            print(f"Attempting to stop container {container_id}")
            container_cleanup_success = container_manager.stop_task_container(task_id, timeout=10)
        
        # 4. 컨테이너 이름 패턴으로도 정리 시도 (백업 방법)
        if not container_cleanup_success:
            print("Attempting container cleanup by name pattern")
            try:
                # task_id의 앞 8자리를 사용하여 컨테이너 이름 패턴 매칭
                task_short_id = task_id[:8]
                container_name_pattern = f"*task-{task_short_id}*"
                container_manager.stop_container_by_name(container_name_pattern)
            except Exception as pattern_error:
                print(f"Pattern-based container cleanup failed: {pattern_error}")
        
        # 5. 작업 상태 확인 및 DB 업데이트
        task = get_task_info(task_id)
        print(f"Task info after revoke: {task}")

        task_status = task.get("status")
        if task_status == 'REVOKED':
            return {
                "message": "Task Revoked Successfully", 
                "task_id": task_id,
                "container_cleanup": container_cleanup_success
            }
        else:
            # 태스크 상태를 'REVOKED'로 강제 업데이트
            crud_task.end_task(current_user.id, task_id, datetime.now(), 'REVOKED')
            print(f"Forced task status update to REVOKED for task {task_id}")
            
            # 업데이트 후 다시 확인
            task = get_task_info(task_id)
            task_status = task.get("status")
            
            if task_status == 'REVOKED':
                return {
                    "message": "Task Revoked Successfully (Forced Update)", 
                    "task_id": task_id,
                    "container_cleanup": container_cleanup_success
                }
            else:
                return {
                    "message": "Task Revoke Completed with Warnings", 
                    "task_id": task_id,
                    "warning": "Task status update may be delayed",
                    "container_cleanup": container_cleanup_success
                }
    
    except Exception as e:
        print(f"Error during task revocation: {str(e)}")
        
        # 오류 발생 시에도 컨테이너 정리 시도
        try:
            container_manager.stop_task_container(task_id, timeout=5)
        except Exception as cleanup_error:
            print(f"Emergency container cleanup failed: {cleanup_error}")
        
        return {
            "message": "Task Revoke Failed", 
            "task_id": task_id,
            "error": str(e)
        }

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

@router.get("/containers/status")
def get_container_status(
    current_user: models.User = Depends(dep.get_current_active_user)
) -> dict:
    """
    Get current container status for debugging
    """
    try:
        import docker
        client = docker.from_env()
        
        # 실행 중인 플러그인 컨테이너 조회
        containers = client.containers.list(
            filters={"label": "container.type=plugin-execution"}
        )
        
        container_info = []
        for container in containers:
            labels = container.labels
            
            # 환경변수에서 CELERY_TASK_ID 찾기
            env_task_id = "unknown"
            env_vars = container.attrs.get('Config', {}).get('Env', [])
            for env in env_vars:
                if env.startswith('CELERY_TASK_ID='):
                    env_task_id = env.split('=', 1)[1]
                    break
            
            container_info.append({
                "id": container.id[:12],
                "name": container.name,
                "status": container.status,
                "task_id_label": labels.get("celery.task_id", "unknown"),
                "task_id_env": env_task_id,
                "plugin_name": labels.get("plugin.name", "unknown"),
                "created": str(container.attrs["Created"]),
                "is_tracked": container.id in container_manager._container_tasks
            })
        
        # 컨테이너 매니저 상태 추가
        manager_status = container_manager.get_status()
        
        return {
            "running_containers": len(container_info),
            "containers": container_info,
            "container_manager": manager_status,
            "orphaned_containers": [
                c for c in container_info 
                if not c["is_tracked"] and c["task_id_label"] != "unknown"
            ]
        }
        
    except Exception as e:
        return {
            "error": f"Failed to get container status: {str(e)}"
        }

@router.get("/logs/{task_id}")
def get_task_logs(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
) -> dict:
    """
    특정 task의 로그 파일들을 조회합니다.
    """
    try:
        # task_id로 task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        
        # 권한 확인 (해당 유저의 task인지 확인)
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")
        
        # 로그 폴더 경로 구성
        logs_folder_path = f"./user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}/logs"
        
        if not os.path.exists(logs_folder_path):
            return {
                "message": "Logs folder not found", 
                "logs": [],
                "task_info": {
                    "task_id": task.task_id,
                    "workflow_id": task.workflow_id,
                    "algorithm_id": task.algorithm_id,
                    "status": task.status
                }
            }
        
        # 로그 파일들 읽기
        log_files = []
        logs_path = Path(logs_folder_path)
        
        for log_file_path in logs_path.glob("*"):
            if log_file_path.is_file():
                try:
                    with open(log_file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    log_files.append({
                        "filename": log_file_path.name,
                        "content": content,
                        "size": log_file_path.stat().st_size,
                        "modified_time": str(log_file_path.stat().st_mtime)
                    })
                except Exception as e:
                    # 파일 읽기 실패 시에도 파일 정보는 포함
                    log_files.append({
                        "filename": log_file_path.name,
                        "content": f"Error reading file: {str(e)}",
                        "size": log_file_path.stat().st_size,
                        "modified_time": str(log_file_path.stat().st_mtime)
                    })
        
        # run.log 파일을 맨 앞으로 정렬
        log_files.sort(key=lambda x: (x["filename"] != "run.log", x["filename"]))
        
        return {
            "message": "Logs retrieved successfully",
            "logs": log_files,
            "task_info": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "status": task.status,
                "start_time": str(task.start_time),
                "end_time": str(task.end_time) if task.end_time else None
            }
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving logs: {str(e)}")
