from fastapi import APIRouter, Depends
from fastapi.responses import StreamingResponse
from sse_starlette.sse import EventSourceResponse
from typing import Any, List, Dict, Optional
from sqlalchemy.orm import Session

from app.task.schemas import TaskMonitoringItem
from app.auth import deps as dep
from app import models
from app.task import service
from app.task import sse

# --- Test-patch compatibility re-exports ---------------------------------
# The business logic now lives in ``app.task.service`` and the SSE streaming in
# ``app.task.sse``; the router only wires HTTP concerns to them. Existing tests
# patch these symbols on the service/sse namespaces (e.g.
# ``app.task.service.parse_snakefile_native``, ``app.task.sse.get_task_info``).
# A couple of symbols are re-exported here for historical patch compatibility.
# noqa: F401 imports below are intentional re-exports for test compatibility.
from app.worker.utils import get_task_info  # noqa: F401

router = APIRouter()


@router.get("/resources")
async def get_resources():
    """리소스 할당 현황 및 실행 중 작업 목록 조회"""
    return service.get_resources()


@router.get("/info/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task.
    SSE 연결은 1시간 타임아웃 및 정리 로직이 적용됩니다.
    """
    return EventSourceResponse(sse.stream_task_status(task_id))


@router.get("/status/{task_id}")
async def get_task_status_simple(task_id: str) -> dict:
    """
    SSE 재연결 판단용 경량 상태 확인 엔드포인트.
    Celery AsyncResult 상태만 반환. 인증 불필요 (SSE 엔드포인트와 동일 패턴).
    """
    return service.get_task_status_simple(task_id=task_id)


@router.get("/monitoring", response_model=List[TaskMonitoringItem])
async def get_task_monitoring(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> List[TaskMonitoringItem]:
    """
    Return the status of all User Tasks with workflow information.
    Excludes plugin build tasks and returns only workflow-related tasks.
    Uses optimized queries to prevent N+1 query problems.
    """
    return service.get_task_monitoring(db=db, current_user=current_user)


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
    return service.revoke_task(db=db, current_user=current_user, task_id=task_id)


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
    return service.delete_task(db=db, current_user=current_user, task_id=task_id)


@router.get("/containers/status")
def get_container_status(
    current_user: models.User = Depends(dep.get_current_active_user)
) -> dict:
    """
    Get current container status for debugging
    """
    return service.get_container_status(current_user=current_user)


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
    return service.get_task_logs(db=db, current_user=current_user, task_id=task_id)


@router.get("/logs/{task_id}/export/json")
def export_task_logs_json(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
) -> StreamingResponse:
    """
    특정 task의 모든 로그 파일들을 JSON 형식으로 export합니다.
    JSON 형태: {"filename1.log": "content1", "filename2.log": "content2"}
    """
    return service.export_task_logs_json(db=db, current_user=current_user, task_id=task_id)


@router.get("/logs/{task_id}/export/txt/{filename}")
def export_task_log_txt(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str,
    filename: str
) -> StreamingResponse:
    """
    특정 task의 개별 로그 파일을 TXT 형식으로 export합니다.
    """
    return service.export_task_log_txt(
        db=db, current_user=current_user, task_id=task_id, filename=filename
    )


@router.get("/{task_id}/dag-structure")
async def get_dag_structure(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
) -> Dict[str, Any]:
    """
    특정 Task의 Snakefile을 파싱하여 DAG 구조를 반환합니다.
    """
    return service.get_dag_structure(db=db, current_user=current_user, task_id=task_id)


@router.get("/{task_id}/rule-status")
async def get_rule_status(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str,
    actual_status: Optional[str] = None
) -> Dict[str, Any]:
    """
    특정 Task의 각 Rule별 실행 상태를 로그 파일 기반으로 추적하여 반환합니다.
    """
    return service.get_rule_status(
        db=db, current_user=current_user, task_id=task_id, actual_status=actual_status
    )


@router.get("/{task_id}/enhanced-progress")
async def get_enhanced_progress(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
) -> Dict[str, Any]:
    """
    특정 Task의 향상된 진행률 정보를 반환합니다.

    향상된 기능:
    - 타이밍 기반 분석
    - 예상 완료 시간
    - 병목 현상 감지
    - 정체 상황 감지
    - 상세 진행률 계산
    """
    return service.get_enhanced_progress(db=db, current_user=current_user, task_id=task_id)


@router.get("/{task_id}/rule-logs/{rule_name}")
async def get_rule_logs(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str,
    rule_name: str
) -> Dict[str, Any]:
    """
    특정 Rule의 로그 파일들(stdout, stderr)을 조회합니다.
    """
    return service.get_rule_logs(
        db=db, current_user=current_user, task_id=task_id, rule_name=rule_name
    )


@router.get("/cache/stats")
def get_dag_cache_stats(
    current_user: models.User = Depends(dep.get_current_active_user)
) -> dict:
    """
    DAG 파싱 캐시 통계 정보 조회 (관리자 전용)
    """
    return service.get_dag_cache_stats(current_user=current_user)


@router.delete("/cache/clear")
def clear_dag_caches(
    current_user: models.User = Depends(dep.get_current_active_user)
) -> dict:
    """
    모든 DAG 파싱 캐시 초기화 (관리자 전용)
    """
    return service.clear_dag_caches(current_user=current_user)


@router.get("/{task_id}/execution-manifest")
def get_execution_manifest(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    task_id: str
) -> StreamingResponse:
    """
    SUCCESS 상태이고 Analysis 타입 플러그인인 Task의 execution manifest를 JSON 형식으로 다운로드합니다.
    분석 재현을 위한 포괄적인 정보를 포함합니다.

    포함 정보:
    - Task, Plugin, Workflow 메타데이터
    - 모든 로그 파일 (run.log, *.stdout, *.stderr 등)
    - Snakefile 내용
    - Plugin metadata.json 내용
    - meta.yml 파일 내용 (있는 경우)
    """
    return service.get_execution_manifest(db=db, current_user=current_user, task_id=task_id)
