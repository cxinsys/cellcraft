from fastapi import APIRouter, Depends, Query
from typing import Any
from sqlalchemy.orm import Session

from app.auth import deps as dep
from app import models
from app.admin import service

router = APIRouter()


@router.get("/users", response_model=Any)
def get_filtered_users(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
):
    return service.get_filtered_users(
        db=db, current_user=current_user, amount=amount, page_num=page_num,
        sort=sort, order=order, searchTerm=searchTerm,
    )


@router.get("/users_count", response_model=Any)
def get_users_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_users_count(db=db, current_user=current_user)


@router.get("/files", response_model=Any)
def get_filtered_files(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):
    return service.get_filtered_files(
        db=db, current_user=current_user, amount=amount, page_num=page_num,
        sort=sort, order=order, searchTerm=searchTerm,
    )


@router.get("/files_count", response_model=Any)
def get_files_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_files_count(db=db, current_user=current_user)


@router.get("/workflows", response_model=Any)
def get_filtered_workflows(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
):
    return service.get_filtered_workflows(
        db=db, current_user=current_user, amount=amount, page_num=page_num,
        sort=sort, order=order, searchTerm=searchTerm,
    )


@router.get("/workflows_count", response_model=Any)
def get_workflows_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_workflows_count(db=db, current_user=current_user)


@router.get("/tasks", response_model=Any)
def get_filtered_tasks(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):
    return service.get_filtered_tasks(
        db=db, current_user=current_user, amount=amount, page_num=page_num,
        sort=sort, order=order, searchTerm=searchTerm,
    )


@router.get("/tasks_count", response_model=Any)
def get_tasks_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_tasks_count(db=db, current_user=current_user)


@router.get("/plugins", response_model=Any)
def get_filtered_plugins(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):
    return service.get_filtered_plugins(
        db=db, current_user=current_user, amount=amount, page_num=page_num,
        sort=sort, order=order, searchTerm=searchTerm,
    )


@router.get("/plugins_count", response_model=Any)
def get_plugins_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_plugins_count(db=db, current_user=current_user)


@router.put("/users/{user_id}", response_model=Any)
def update_user(
    user_id: int,
    user_data: dict,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.update_user(
        db=db, current_user=current_user, user_id=user_id, user_data=user_data
    )


@router.delete("/users/{user_id}")
def delete_user(
    user_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.delete_user(db=db, current_user=current_user, user_id=user_id)


@router.delete("/files/{file_id}")
def delete_file(
    file_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.delete_file(db=db, current_user=current_user, file_id=file_id)


@router.delete("/workflows/{workflow_id}")
def delete_workflow(
    workflow_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.delete_workflow(db=db, current_user=current_user, workflow_id=workflow_id)


@router.post("/tasks/{task_id}/cancel")
def cancel_task(
    task_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.cancel_task(db=db, current_user=current_user, task_id=task_id)


@router.post("/plugins/{plugin_id}/install-dependencies")
def install_plugin_dependencies(
    plugin_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.install_plugin_dependencies(
        db=db, current_user=current_user, plugin_id=plugin_id
    )


# =============================================================================
# Plugin Synchronization Management API Endpoints
# =============================================================================

@router.get("/plugins/sync/status")
def get_plugin_sync_status(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """Get current plugin synchronization status"""
    return service.get_plugin_sync_status(current_user=current_user)


@router.post("/plugins/sync")
def sync_plugins_from_repository(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """Manually trigger plugin synchronization from repository to database"""
    return service.sync_plugins_from_repository(current_user=current_user)


@router.get("/plugins/consistency")
def get_plugin_consistency_report(
    current_user: models.User = Depends(dep.get_current_active_user),
    format: str = Query("json", description="Report format: json or text")
):
    """Get plugin version consistency report"""
    return service.get_plugin_consistency_report(current_user=current_user, format=format)


@router.get("/plugins/branches")
def get_available_plugin_branches(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """Get list of available plugin repository branches"""
    return service.get_available_plugin_branches(current_user=current_user)


@router.post("/plugins/branch/{branch}")
def switch_plugin_branch(
    branch: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """Switch plugin repository branch and synchronize database"""
    return service.switch_plugin_branch(current_user=current_user, branch=branch)


@router.get("/plugins/consistency/quick")
def quick_consistency_check(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """Quick consistency check (returns boolean only)"""
    return service.quick_consistency_check(current_user=current_user)


# =============================================================================
# Container Manager API Endpoints (메모리 누수 방지)
# =============================================================================

@router.get("/container-manager/status")
def get_container_manager_status(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    Container Manager 상태 조회

    Returns:
        - tracked_tasks: 추적 중인 작업 수
        - cleanup_in_progress_count: 정리 중인 컨테이너 수
        - task_container_mapping: 작업-컨테이너 매핑
        - container_task_mapping: 컨테이너-작업 매핑
        - cleanup_in_progress: 정리 중인 컨테이너 ID 목록
    """
    return service.get_container_manager_status(current_user=current_user)


@router.post("/container-manager/cleanup")
def cleanup_container_manager(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    Container Manager stale 매핑 수동 정리

    실제로 존재하지 않는 Docker 컨테이너에 대한 매핑을 정리합니다.
    이 작업은 장시간 운영 후 메모리 누수 방지를 위해 주기적으로 수행할 수 있습니다.

    Returns:
        - task_containers: 정리된 작업-컨테이너 매핑 수
        - cleanup_set: 정리된 cleanup_in_progress 항목 수
        - errors: 발생한 오류 목록 (있는 경우)
    """
    return service.cleanup_container_manager(current_user=current_user)


@router.post("/container-manager/force-clear")
def force_clear_container_manager(
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    Container Manager 모든 매핑 강제 정리 (비상용)

    WARNING: 이 작업은 모든 매핑을 즉시 삭제합니다.
    실행 중인 작업의 추적이 불가능해질 수 있으므로,
    시스템 재시작 전이나 비상 상황에서만 사용하세요.

    Returns:
        - cleared_tasks: 삭제된 작업-컨테이너 매핑 수
        - cleared_cleanup_markers: 삭제된 cleanup 마커 수
    """
    return service.force_clear_container_manager(current_user=current_user)
