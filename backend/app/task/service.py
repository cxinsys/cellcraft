"""
Task domain service layer.

Extracted from ``task/router.py`` in PR-7 (Phase 3c). This module holds the
business logic that previously lived inline in the endpoints: resource status,
lightweight task status, task monitoring, revoke/delete, container status, log
retrieval + export, DAG structure / rule-status / enhanced-progress / rule-logs
parsing, DAG cache stats/clear, and execution manifest generation. The router
now only parses requests, resolves dependencies, calls a service function, and
returns its result. SSE streaming lives in ``task/sse.py``.

Transitional exception policy (PR-7): the global domain-exception layer arrives
in PR-8. Until then, service functions may ``raise HTTPException`` directly
(copy-move fidelity first, exception translation later). Response shapes — JSON
keys, status codes, and error message strings — are preserved exactly.

Snakefile parser modules (``workflow/compiler/dag_parser``,
``workflow/compiler/native_parser``) are invoked, never modified (PR-10).
"""
import io
import json
import logging
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from fastapi.responses import StreamingResponse
from sqlalchemy.orm import Session

from app import models
from app.core.enums import PluginType
from app.core.exceptions import ValidationFailedError, ExternalServiceError
from app.shared.docker import container_manager
from app.task import crud as crud_task
from app.task.schemas import TaskMonitoringItem, PluginInfo
from app.workflow import crud as crud_workflow
from app.workflow.compiler.dag_parser import (
    SnakemakeDAGParser, SnakemakeRuleStatusTracker
)
from app.workflow.compiler.native_parser import (
    parse_snakefile_native, SnakemakeNativeError
)
from app.worker.utils import get_task_info

logger = logging.getLogger(__name__)


def get_resources() -> dict:
    """리소스 할당 현황 및 실행 중 작업 목록 조회."""
    from app.shared.resources import get_resource_status

    status = get_resource_status()
    if status is None:
        raise ExternalServiceError("Resource monitoring unavailable", status_code=503)
    return status


def get_task_status_simple(*, task_id: str) -> dict:
    """SSE 재연결 판단용 경량 상태 확인 (Celery AsyncResult 상태만 반환)."""
    if not task_id or not task_id.strip():
        raise ValidationFailedError("Invalid task_id")

    task = get_task_info(task_id)
    return {
        "task_id": task_id,
        "status": task.get("task_status", "UNKNOWN")
    }


def get_task_monitoring(*, db: Session, current_user: models.User) -> List[TaskMonitoringItem]:
    """Return all of a user's workflow tasks (plugin_build filtered out).

    Uses eager-loaded plugin/workflow data (single query) to avoid N+1.
    """
    # Get user tasks with eager-loaded plugin and workflow data (single query)
    user_tasks = crud_task.get_user_task_with_plugin(db, current_user.id)

    # Return empty list with 200 OK if no tasks (RESTful response)
    if not user_tasks:
        return []

    # Filter out plugin build tasks, keep only workflow-related tasks
    workflow_tasks = [task for task in user_tasks if task.task_type != "plugin_build"]

    # Build response using eager-loaded data (no additional queries needed)
    tasks_with_workflow = []
    for task in workflow_tasks:
        # task_title generation logic
        task_title = None
        if task.plugin_name:
            if task.task_type == 'compile':
                task_title = task.plugin_name
            elif task.task_type == 'visualization':
                task_title = f"{task.plugin_name}-visualization"
            else:
                task_title = task.plugin_name  # default value

        # Plugin information - use eager-loaded data
        plugin_info = None
        if task.plugin:
            plugin_info = PluginInfo(
                version=task.plugin.version,
                source=task.plugin.source,
                plugin_type=task.plugin.plugin_type.value if task.plugin.plugin_type else None
            )

        # Use eager-loaded workflow data (no additional query)
        workflow_title = task.workflows.title if task.workflows else None

        task_item = TaskMonitoringItem(
            id=task.id,
            task_id=task.task_id,
            workflow_id=task.workflow_id,
            user_id=task.user_id,
            status=task.status,
            start_time=task.start_time,
            end_time=task.end_time,
            workflow_title=workflow_title,
            task_title=task_title,
            plugin_name=task.plugin_name,
            task_type=task.task_type,
            plugin=plugin_info
        )
        tasks_with_workflow.append(task_item)

    return tasks_with_workflow


def revoke_task(*, db: Session, current_user: models.User, task_id: str) -> dict:
    """Revoke a Celery task and clean up its containers + results folder."""
    try:
        # 지연 import 유지 (사유): 테스트 patch 호환(`app.task.service.get_celery_app`가
        # 아닌 호출부 재바인딩 방식) + worker가 task.service를 import 해도 web 계층
        # (app.main)이 로드되지 않도록 하기 위함 (import-linter no-web-in-worker).
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

        # 5. Results 폴더 정리 (재실행 시 Snakemake 스킵 방지)
        results_cleanup_success = False
        try:
            task_record = crud_task.get_task_by_task_id(db, task_id)
            if task_record:
                if task_record.task_type == 'visualization':
                    task_folder = f"./user/{current_user.username}/workflow_{task_record.workflow_id}/visualization_{task_record.algorithm_id}"
                else:
                    task_folder = f"./user/{current_user.username}/workflow_{task_record.workflow_id}/algorithm_{task_record.algorithm_id}"

                from app.shared.archive import cleanup_task_results
                cleanup_result = cleanup_task_results(Path(task_folder), preserve_folder=True)
                results_cleanup_success = cleanup_result["success"]

                if results_cleanup_success:
                    print(f"Results cleanup completed: {len(cleanup_result.get('files_removed', []))} files removed")
                else:
                    print(f"Results cleanup failed: {cleanup_result.get('error')}")
        except Exception as cleanup_error:
            print(f"Results cleanup failed: {cleanup_error}")

        # 6. 작업 상태 확인 및 DB 업데이트
        task = get_task_info(task_id)
        print(f"Task info after revoke: {task}")

        task_status = task.get("status")
        if task_status == 'REVOKED':
            return {
                "message": "Task Revoked Successfully",
                "task_id": task_id,
                "container_cleanup": container_cleanup_success,
                "results_cleanup": results_cleanup_success
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
                    "container_cleanup": container_cleanup_success,
                    "results_cleanup": results_cleanup_success
                }
            else:
                return {
                    "message": "Task Revoke Completed with Warnings",
                    "task_id": task_id,
                    "warning": "Task status update may be delayed",
                    "container_cleanup": container_cleanup_success,
                    "results_cleanup": results_cleanup_success
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


def delete_task(*, db: Session, current_user: models.User, task_id: str) -> dict:
    """Delete a user's task record."""
    task = crud_task.delete_user_task(db, current_user.id, task_id)
    return {"message": "Task Deleted", "task_id": task_id}


def get_container_status(*, current_user: models.User) -> dict:
    """Get current plugin-execution container status for debugging."""
    client = None  # Docker 클라이언트 누수 방지를 위해 초기화
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
    finally:
        # Docker 클라이언트 연결 해제 - TCP 소켓 누수 방지
        if client:
            try:
                client.close()
            except Exception:
                pass


def get_task_logs(*, db: Session, current_user: models.User, task_id: str) -> dict:
    """특정 task의 로그 파일들을 조회합니다."""
    from app.file.storage import construct_logs_path, is_safe_path, get_logs_base_path

    try:
        # task_id로 task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인 (해당 유저의 task인지 확인)
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")

        # 로그 폴더 경로 구성 (centralized utility 사용)
        # task_id를 전달하여 아카이브된 로그가 있으면 해당 경로 사용
        logs_folder_path = construct_logs_path(
            username=current_user.username,
            workflow_id=task.workflow_id,
            algorithm_id=task.algorithm_id,
            task_type=task.task_type,
            task_id=task.task_id
        )

        # 보안 검증: path traversal 방지
        if not is_safe_path(logs_folder_path, get_logs_base_path()):
            raise HTTPException(status_code=400, detail="Invalid path")

        if not os.path.exists(logs_folder_path):
            raise HTTPException(status_code=404, detail="Logs folder not found")

        # 로그 파일들 읽기
        log_files = []
        logs_path = Path(logs_folder_path)

        for log_file_path in logs_path.glob("*"):
            if log_file_path.is_file():
                # 각 파일도 안전성 검증
                if not is_safe_path(str(log_file_path), get_logs_base_path()):
                    continue  # Skip unsafe files (e.g., symlinks outside base dir)

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


def export_task_logs_json(*, db: Session, current_user: models.User, task_id: str) -> StreamingResponse:
    """특정 task의 모든 로그 파일들을 JSON 형식으로 export합니다."""
    from app.file.storage import construct_logs_path, is_safe_path, get_logs_base_path

    try:
        # task_id로 task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인 (해당 유저의 task인지 확인)
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")

        # 로그 폴더 경로 구성
        # task_id를 전달하여 아카이브된 로그가 있으면 해당 경로 사용
        logs_folder_path = construct_logs_path(
            username=current_user.username,
            workflow_id=task.workflow_id,
            algorithm_id=task.algorithm_id,
            task_type=task.task_type,
            task_id=task.task_id
        )

        # 보안 검증
        if not is_safe_path(logs_folder_path, get_logs_base_path()):
            raise HTTPException(status_code=400, detail="Invalid path")

        if not os.path.exists(logs_folder_path):
            raise HTTPException(status_code=404, detail="Logs folder not found")

        # 모든 로그 파일을 JSON 형태로 수집
        logs_data = {}
        logs_path = Path(logs_folder_path)

        for log_file_path in logs_path.glob("*"):
            if log_file_path.is_file():
                # 보안 검증: 각 파일이 안전한지 확인
                if not is_safe_path(str(log_file_path), get_logs_base_path()):
                    continue  # Skip unsafe files

                try:
                    with open(log_file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    logs_data[log_file_path.name] = content
                except Exception as e:
                    # 파일 읽기 실패 시 에러 정보 포함
                    logs_data[log_file_path.name] = f"Error reading file: {str(e)}"

        if not logs_data:
            raise HTTPException(status_code=404, detail="No log files found")

        # JSON으로 변환
        json_content = json.dumps(logs_data, indent=2, ensure_ascii=False)

        # 파일명 생성 (타임스탬프 포함)
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"task_{task_id[:8]}_logs_{timestamp}.json"

        # StreamingResponse로 다운로드 제공
        return StreamingResponse(
            io.StringIO(json_content),
            media_type="application/json",
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error exporting logs as JSON: {str(e)}")


def export_task_log_txt(*, db: Session, current_user: models.User, task_id: str,
                        filename: str) -> StreamingResponse:
    """특정 task의 개별 로그 파일을 TXT 형식으로 export합니다."""
    from app.file.storage import (
        construct_logs_path,
        is_safe_path,
        get_logs_base_path,
        validate_export_filename
    )

    try:
        # 파일명 보안 검증 (path traversal 방지)
        if not validate_export_filename(filename, task_id):
            raise HTTPException(status_code=400, detail="Invalid filename")

        # task_id로 task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인 (해당 유저의 task인지 확인)
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")

        # 로그 폴더 경로 구성
        # task_id를 전달하여 아카이브된 로그가 있으면 해당 경로 사용
        logs_folder_path = construct_logs_path(
            username=current_user.username,
            workflow_id=task.workflow_id,
            algorithm_id=task.algorithm_id,
            task_type=task.task_type,
            task_id=task.task_id
        )

        # 보안 검증: 폴더 경로
        if not is_safe_path(logs_folder_path, get_logs_base_path()):
            raise HTTPException(status_code=400, detail="Invalid path")

        # 요청된 로그 파일 경로
        log_file_path = Path(logs_folder_path) / filename

        # 보안 검증: 최종 파일 경로
        if not is_safe_path(str(log_file_path), get_logs_base_path()):
            raise HTTPException(status_code=400, detail="Invalid file path")

        if not log_file_path.exists():
            raise HTTPException(status_code=404, detail=f"Log file '{filename}' not found")

        if not log_file_path.is_file():
            raise HTTPException(status_code=400, detail=f"'{filename}' is not a valid file")

        # 파일 읽기
        try:
            with open(log_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Error reading file: {str(e)}")

        # 다운로드용 파일명 생성 (타임스탬프 포함)
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = Path(filename).stem
        extension = Path(filename).suffix or ".txt"
        download_filename = f"task_{task_id[:8]}_{base_name}_{timestamp}{extension}"

        # StreamingResponse로 다운로드 제공
        return StreamingResponse(
            io.StringIO(content),
            media_type="text/plain",
            headers={"Content-Disposition": f"attachment; filename={download_filename}"}
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error exporting log file as TXT: {str(e)}")


def get_dag_structure(*, db: Session, current_user: models.User, task_id: str) -> Dict[str, Any]:
    """특정 Task의 Snakefile을 파싱하여 DAG 구조를 반환합니다."""
    try:
        # 입력 검증
        if not task_id or not task_id.strip():
            raise HTTPException(status_code=400, detail="Invalid task_id")

        # Task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            logger.warning(f"Task not found: {task_id}")
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인
        if task.user_id != current_user.id:
            logger.warning(f"Access denied for user {current_user.id} to task {task_id}")
            raise HTTPException(status_code=403, detail="Access denied")

        # Snakefile 경로 구성
        try:
            if task.task_type == 'visualization':
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}/Snakefile"
            else:
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}/Snakefile"
        except AttributeError as e:
            logger.error(f"Missing task attributes for DAG structure: {e}")
            raise HTTPException(status_code=400, detail="Invalid task data")

        logger.info(f"Parsing DAG structure for task {task_id}, snakefile: {snakefile_path}")

        # Snakefile 파싱 (네이티브 파서 우선 사용, 실패 시 레거시 파서로 fallback)
        dag_data = None
        parsing_method = "unknown"

        try:
            # 1단계: Snakemake 네이티브 파서 시도 (실시간 파싱)
            logger.info(f"Attempting native Snakemake parsing for task {task_id}")
            dag_data = parse_snakefile_native(snakefile_path, method='auto')
            parsing_method = f"native_{dag_data.get('method', 'unknown')}"

            # DAG 데이터 유효성 검증
            if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                raise ValueError("Invalid DAG data structure from native parser")

            logger.info(f"Successfully parsed {len(dag_data.get('nodes', []))} rules using native parser ({parsing_method}) for task {task_id}")

        except (SnakemakeNativeError, ValueError, Exception) as native_error:
            logger.warning(f"Native parser failed for task {task_id}: {native_error}")

            try:
                # 2단계: 레거시 파서로 fallback
                logger.info(f"Falling back to legacy parser for task {task_id}")
                parser = SnakemakeDAGParser()
                dag_data = parser.parse_snakefile_with_logs(snakefile_path)
                parsing_method = "legacy"

                # DAG 데이터 유효성 검증
                if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                    logger.error(f"Invalid DAG data structure from legacy parser for task {task_id}")
                    raise HTTPException(status_code=500, detail="Invalid DAG structure")

                logger.info(f"Successfully parsed {len(dag_data.get('nodes', []))} rules using legacy parser for task {task_id}")

            except Exception as legacy_error:
                # 두 파서 모두 실패한 경우
                logger.error(f"Both native and legacy parsers failed for task {task_id}. Native: {native_error}, Legacy: {legacy_error}")
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to parse workflow with both native and legacy parsers. Last error: {str(legacy_error)}"
                )

        # 파싱 방법 정보 추가
        if dag_data and isinstance(dag_data, dict):
            dag_data['parsing_method'] = parsing_method

        return {
            "task_info": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "task_type": task.task_type,
                "plugin_name": task.plugin_name,
                "status": task.status
            },
            "dag_structure": dag_data,
            "snakefile_path": snakefile_path
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Unexpected error getting DAG structure for task {task_id}: {e}")
        raise HTTPException(status_code=500, detail="Internal server error while processing DAG structure")


def get_rule_status(*, db: Session, current_user: models.User, task_id: str,
                    actual_status: Optional[str] = None) -> Dict[str, Any]:
    """특정 Task의 각 Rule별 실행 상태를 로그 파일 기반으로 추적하여 반환합니다."""
    try:
        # 입력 검증
        if not task_id or not task_id.strip():
            raise HTTPException(status_code=400, detail="Invalid task_id")

        # Task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            logger.warning(f"Task not found for rule status: {task_id}")
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인
        if task.user_id != current_user.id:
            logger.warning(f"Access denied for user {current_user.id} to task rule status {task_id}")
            raise HTTPException(status_code=403, detail="Access denied")

        # Snakefile 경로 구성
        try:
            if task.task_type == 'visualization':
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}/Snakefile"
            else:
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}/Snakefile"
        except AttributeError as e:
            logger.error(f"Missing task attributes for rule status: {e}")
            raise HTTPException(status_code=400, detail="Invalid task data")

        logger.info(f"Getting rule status for task {task_id}, snakefile: {snakefile_path}")

        # Snakefile 파싱 (네이티브 파서 우선 사용)
        dag_data = None
        parsing_method = "unknown"

        try:
            # 1단계: Snakemake 네이티브 파서 시도 (실시간 파싱)
            dag_data = parse_snakefile_native(snakefile_path, method='auto')
            parsing_method = f"native_{dag_data.get('method', 'unknown')}"

            # DAG 데이터 유효성 검증
            if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                raise ValueError("Invalid DAG data structure from native parser")

        except (SnakemakeNativeError, ValueError, Exception) as native_error:
            logger.warning(f"Native parser failed for rule status {task_id}: {native_error}")

            try:
                # 2단계: 레거시 파서로 fallback
                parser = SnakemakeDAGParser()
                dag_data = parser.parse_snakefile_with_logs(snakefile_path)
                parsing_method = "legacy"

                # DAG 데이터 유효성 검증
                if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                    logger.error(f"Invalid DAG data structure for rule status {task_id}")
                    raise HTTPException(status_code=500, detail="Invalid workflow structure")

            except Exception as legacy_error:
                logger.error(f"Both parsers failed for rule status {task_id}. Native: {native_error}, Legacy: {legacy_error}")
                raise HTTPException(
                    status_code=500,
                    detail=f"Cannot parse workflow: {str(legacy_error)}"
                )

        # 상태 추적기 초기화 및 향상된 진행률 계산
        try:
            # Use actual_status if provided, otherwise use task.status
            effective_task_status = actual_status if actual_status else task.status
            logger.info(f"Using task status for rule tracking: {effective_task_status} (actual_status: {actual_status}, db_status: {task.status})")

            tracker = SnakemakeRuleStatusTracker(
                workflow_path=snakefile_path,
                task_status=effective_task_status
            )
            tracker.rules = dag_data['nodes']

            # 향상된 진행률 정보 가져오기
            enhanced_progress = tracker.get_enhanced_progress_info()

            # 기본 진행률 정보 추출
            basic_progress = enhanced_progress.get('basic_progress', {})
            rule_statuses = enhanced_progress.get('rule_statuses', {})

            # 기존 변수들 설정 (하위 호환성)
            total_rules = basic_progress.get('total_rules', len(dag_data.get('execution_sequence', [])))
            completed_rules = basic_progress.get('completed_rules', 0)
            failed_rules = basic_progress.get('failed_rules', 0)
            running_rules = basic_progress.get('running_rules', 0)
            pending_rules = basic_progress.get('pending_rules', 0)
            progress_percentage = basic_progress.get('percentage', 0.0)

            logger.info(f"Enhanced rule status calculated for task {task_id}: {completed_rules}/{total_rules} completed")
            logger.info(f"Rule statuses for task {task_id}: {rule_statuses}")
            logger.info(f"Rule statuses type: {type(rule_statuses)}, keys: {list(rule_statuses.keys()) if isinstance(rule_statuses, dict) else 'Not a dict'}")

        except Exception as e:
            logger.error(f"Error tracking rule statuses for task {task_id}: {e}")
            # 기본 상태로 모든 룰을 pending으로 설정
            rule_statuses = {rule['id']: 'pending' for rule in dag_data.get('nodes', [])}

            # 폴백 진행률 계산
            try:
                total_rules = len(dag_data.get('execution_sequence', []))
                completed_rules = failed_rules = running_rules = 0
                pending_rules = total_rules
                progress_percentage = 0
                enhanced_progress = None
            except Exception:
                total_rules = len(rule_statuses)
                completed_rules = failed_rules = running_rules = 0
                pending_rules = total_rules
                progress_percentage = 0
                enhanced_progress = None

        # 응답 구성
        response = {
            "task_info": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "task_type": task.task_type,
                "plugin_name": task.plugin_name,
                "status": task.status
            },
            "rule_statuses": rule_statuses,
            "execution_sequence": dag_data.get('execution_sequence', []),
            "progress": {
                "total_rules": total_rules,
                "completed_rules": completed_rules,
                "failed_rules": failed_rules,
                "running_rules": running_rules,
                "pending_rules": pending_rules,
                "percentage": round(progress_percentage, 1)
            },
            "snakefile_path": snakefile_path
        }

        logger.info(f"Final API response rule_statuses for task {task_id}: {response['rule_statuses']}")
        logger.info(f"Final API response type check - rule_statuses is dict: {isinstance(response['rule_statuses'], dict)}")

        # 향상된 진행률 정보 추가 (사용 가능한 경우)
        if 'enhanced_progress' in locals() and enhanced_progress:
            response["enhanced_progress"] = {
                "timing_info": enhanced_progress.get('timing_info', {}),
                "estimated_completion": enhanced_progress.get('estimated_completion', {}),
                "bottleneck_analysis": enhanced_progress.get('bottleneck_analysis', {}),
                "completion_percentage": enhanced_progress.get('basic_progress', {}).get('completion_percentage', 0.0),
                "is_stalled": enhanced_progress.get('basic_progress', {}).get('is_stalled', False)
            }

        return response

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error getting rule status: {str(e)}")


def get_enhanced_progress(*, db: Session, current_user: models.User, task_id: str) -> Dict[str, Any]:
    """특정 Task의 향상된 진행률 정보를 반환합니다."""
    try:
        # 입력 검증
        if not task_id or not task_id.strip():
            raise HTTPException(status_code=400, detail="Invalid task_id")

        # Task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            logger.warning(f"Task not found for enhanced progress: {task_id}")
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인
        if task.user_id != current_user.id:
            logger.warning(f"Access denied for user {current_user.id} to enhanced progress {task_id}")
            raise HTTPException(status_code=403, detail="Access denied")

        # Snakefile 경로 구성
        try:
            if task.task_type == 'visualization':
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}/Snakefile"
            else:
                snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}/Snakefile"
        except AttributeError as e:
            logger.error(f"Missing task attributes for enhanced progress: {e}")
            raise HTTPException(status_code=400, detail="Invalid task data")

        logger.info(f"Getting enhanced progress for task {task_id}, snakefile: {snakefile_path}")

        # Snakefile 파싱 및 향상된 진행률 계산 (네이티브 파서 우선 사용)
        dag_data = None

        try:
            # 1단계: Snakemake 네이티브 파서 시도 (실시간 파싱)
            dag_data = parse_snakefile_native(snakefile_path, method='auto')

            # DAG 데이터 유효성 검증
            if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                raise ValueError("Invalid DAG data structure from native parser")

        except (SnakemakeNativeError, ValueError, Exception) as native_error:
            logger.warning(f"Native parser failed for enhanced progress {task_id}: {native_error}")

            try:
                # 2단계: 레거시 파서로 fallback
                parser = SnakemakeDAGParser()
                dag_data = parser.parse_snakefile_with_logs(snakefile_path)

                # DAG 데이터 유효성 검증
                if not isinstance(dag_data, dict) or 'nodes' not in dag_data:
                    logger.error(f"Invalid DAG data structure for enhanced progress {task_id}")
                    raise HTTPException(status_code=500, detail="Invalid workflow structure")

            except Exception as legacy_error:
                logger.error(f"Both parsers failed for enhanced progress {task_id}. Native: {native_error}, Legacy: {legacy_error}")
                raise HTTPException(
                    status_code=500,
                    detail=f"Cannot parse workflow: {str(legacy_error)}"
                )

        # 상태 추적기 초기화
        tracker = SnakemakeRuleStatusTracker(
            workflow_path=snakefile_path,
            task_status=task.status
        )
        tracker.rules = dag_data['nodes']

        # 향상된 진행률 정보 계산
        enhanced_progress = tracker.get_enhanced_progress_info()

        logger.info(f"Enhanced progress calculated for task {task_id}")

        return {
            "task_info": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "task_type": task.task_type,
                "plugin_name": task.plugin_name,
                "status": task.status
            },
            "enhanced_progress": enhanced_progress,
            "snakefile_path": snakefile_path,
            "timestamp": time.time()
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Unexpected error getting enhanced progress for task {task_id}: {e}")
        raise HTTPException(status_code=500, detail="Internal server error while processing enhanced progress")


def get_rule_logs(*, db: Session, current_user: models.User, task_id: str,
                  rule_name: str) -> Dict[str, Any]:
    """특정 Rule의 로그 파일들(stdout, stderr)을 조회합니다."""
    try:
        # Task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")

        # Snakefile 파싱하여 해당 Rule의 로그 경로 찾기
        if task.task_type == 'visualization':
            snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}/Snakefile"
        else:
            snakefile_path = f"user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}/Snakefile"

        if not os.path.exists(snakefile_path):
            raise HTTPException(status_code=404, detail="Snakefile not found")

        # 네이티브 파서 우선 사용
        try:
            dag_data = parse_snakefile_native(snakefile_path, method='auto')
        except (SnakemakeNativeError, Exception) as native_error:
            logger.warning(f"Native parser failed for rule logs {task_id}: {native_error}")
            # 레거시 파서로 fallback
            parser = SnakemakeDAGParser()
            dag_data = parser.parse_snakefile_with_logs(snakefile_path)

        # 해당 Rule 찾기
        target_rule = None
        for node in dag_data['nodes']:
            if node['id'] == rule_name:
                target_rule = node
                break

        if not target_rule:
            raise HTTPException(status_code=404, detail=f"Rule '{rule_name}' not found")

        # 아카이브된 로그 경로 확인을 위한 기본 경로 구성
        if task.task_type == 'visualization':
            base_task_folder = f"user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}"
        else:
            base_task_folder = f"user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}"

        archived_logs_folder = os.path.join(base_task_folder, "executions", task.task_id, "logs")
        use_archived_logs = os.path.exists(archived_logs_folder)

        # 로그 파일들 읽기
        rule_logs = {}
        for log_type, log_path in target_rule['log_paths'].items():
            # 아카이브된 로그가 있으면 해당 경로로 변환
            actual_log_path = log_path
            if use_archived_logs:
                log_filename = os.path.basename(log_path)
                archived_log_path = os.path.join(archived_logs_folder, log_filename)
                if os.path.exists(archived_log_path):
                    actual_log_path = archived_log_path

            if os.path.exists(actual_log_path):
                try:
                    with open(actual_log_path, 'r', encoding='utf-8') as f:
                        content = f.read()

                    # 파일 정보 포함
                    file_stat = os.stat(actual_log_path)
                    rule_logs[log_type] = {
                        "content": content,
                        "file_path": actual_log_path,
                        "size": file_stat.st_size,
                        "modified_time": file_stat.st_mtime,
                        "exists": True
                    }
                except Exception as e:
                    rule_logs[log_type] = {
                        "content": f"Error reading file: {str(e)}",
                        "file_path": actual_log_path,
                        "size": 0,
                        "modified_time": 0,
                        "exists": True,
                        "error": str(e)
                    }
            else:
                rule_logs[log_type] = {
                    "content": "",
                    "file_path": actual_log_path,
                    "size": 0,
                    "modified_time": 0,
                    "exists": False
                }

        return {
            "task_info": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "task_type": task.task_type
            },
            "rule_info": {
                "rule_name": rule_name,
                "label": target_rule['label'],
                "description": target_rule['description'],
                "inputs": target_rule['inputs'],
                "outputs": target_rule['outputs'],
                "params": target_rule['params']
            },
            "logs": rule_logs
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error getting rule logs: {str(e)}")


def get_dag_cache_stats(*, current_user: models.User) -> dict:
    """DAG 파싱 캐시 통계 정보 조회 (관리자 전용)."""
    try:
        # 레거시 파서의 캐시 통계만 수집 (네이티브 파서는 캐시 제거됨)
        from app.workflow.compiler.dag_parser import get_cache_stats as get_legacy_stats

        legacy_stats = get_legacy_stats()

        return {
            "native_parser_cache": {"message": "Cache removed for simplification"},
            "legacy_parser_cache": legacy_stats,
            "timestamp": time.time(),
            "message": "Cache statistics retrieved successfully"
        }

    except Exception as e:
        logger.error(f"Error getting cache stats: {e}")
        return {
            "error": f"Failed to get cache statistics: {str(e)}",
            "timestamp": time.time()
        }


def clear_dag_caches(*, current_user: models.User) -> dict:
    """모든 DAG 파싱 캐시 초기화 (관리자 전용)."""
    try:
        # 레거시 파서의 캐시만 초기화 (네이티브 파서는 캐시 제거됨)
        from app.workflow.compiler.dag_parser import clear_all_caches

        # 레거시 파서 캐시 초기화
        clear_all_caches()

        logger.info(f"Legacy DAG cache cleared by user {current_user.id} (native cache removed)")

        return {
            "message": "Legacy DAG parsing cache cleared successfully (native cache removed for simplification)",
            "timestamp": time.time(),
            "cleared_by": current_user.username
        }

    except Exception as e:
        logger.error(f"Error clearing caches: {e}")
        return {
            "error": f"Failed to clear caches: {str(e)}",
            "timestamp": time.time()
        }


def get_execution_manifest(*, db: Session, current_user: models.User, task_id: str) -> StreamingResponse:
    """SUCCESS 상태이고 Analysis 타입 플러그인인 Task의 execution manifest를 JSON으로 다운로드합니다."""
    try:
        # task_id로 task 정보 조회
        task = crud_task.get_task_by_task_id(db, task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        # 권한 확인 (해당 유저의 task인지 확인)
        if task.user_id != current_user.id:
            raise HTTPException(status_code=403, detail="Access denied")

        # Task 상태 검증 - SUCCESS 상태만 허용
        if task.status != 'SUCCESS':
            raise HTTPException(
                status_code=400,
                detail=f"Execution manifest is only available for tasks with SUCCESS status. Current status: {task.status}"
            )

        # Plugin 정보 확인 및 타입 검증 - Analysis 타입만 허용
        if not task.plugin:
            raise HTTPException(status_code=400, detail="Plugin information not found for this task")

        if task.plugin.plugin_type != PluginType.ANALYSIS:
            raise HTTPException(
                status_code=400,
                detail=f"Execution manifest is only available for Analysis type plugins. Current plugin type: {task.plugin.plugin_type.value if task.plugin.plugin_type else 'unknown'}"
            )

        # Workflow 정보 조회
        workflow = crud_workflow.get_workflow_by_id(db, task.workflow_id)

        # 로그 폴더 경로 구성 - task_type에 따라 다른 경로 사용
        if task.task_type == 'visualization':
            task_folder_path = f"./user/{current_user.username}/workflow_{task.workflow_id}/visualization_{task.algorithm_id}"
        else:
            task_folder_path = f"./user/{current_user.username}/workflow_{task.workflow_id}/algorithm_{task.algorithm_id}"

        # 아카이브된 로그가 있으면 해당 경로 사용 (task_id별 분리된 로그)
        archived_logs_folder_path = os.path.join(task_folder_path, "executions", task.task_id, "logs")
        if os.path.exists(archived_logs_folder_path):
            logs_folder_path = archived_logs_folder_path
        else:
            logs_folder_path = os.path.join(task_folder_path, "logs")

        # Execution manifest 데이터 구성
        manifest_data = {
            "manifest_info": {
                "format_version": "1.0",
                "generated_at": datetime.now().isoformat(),
                "generated_by": current_user.username,
                "description": "CellCraft execution manifest for analysis reproducibility"
            },
            "task_metadata": {
                "task_id": task.task_id,
                "workflow_id": task.workflow_id,
                "algorithm_id": task.algorithm_id,
                "plugin_name": task.plugin_name,
                "task_type": task.task_type,
                "status": task.status,
                "start_time": str(task.start_time) if task.start_time else None,
                "end_time": str(task.end_time) if task.end_time else None,
                "plugin_image_uri": task.plugin_image_uri
            },
            "plugin_metadata": {
                "id": task.plugin.id,
                "name": task.plugin.name,
                "description": task.plugin.description,
                "author": task.plugin.author,
                "version": task.plugin.version,
                "plugin_type": task.plugin.plugin_type.value if task.plugin.plugin_type else None,
                "source": task.plugin.source,
                "use_gpu": task.plugin.use_gpu,
                "created_at": str(task.plugin.created_at),
                "updated_at": str(task.plugin.updated_at),
                "dependencies": task.plugin.dependencies,
                "drawflow": task.plugin.drawflow,
                "rules": task.plugin.rules
            },
            "workflow_metadata": {
                "id": workflow.id if workflow else None,
                "title": workflow.title if workflow else None,
                "workflow_info": workflow.workflow_info if workflow else None,
                "created_at": str(workflow.created_at) if workflow else None,
                "updated_at": str(workflow.updated_at) if workflow else None
            },
            "execution_files": {
                "logs": {},
                "snakefile": None,
                "plugin_metadata": None,
                "meta_yml": None,
                "results": {}
            }
        }

        # 로그 파일들 수집
        if os.path.exists(logs_folder_path):
            logs_path = Path(logs_folder_path)

            for log_file_path in logs_path.glob("*"):
                if log_file_path.is_file():
                    try:
                        with open(log_file_path, 'r', encoding='utf-8') as f:
                            content = f.read()
                        manifest_data["execution_files"]["logs"][log_file_path.name] = {
                            "content": content,
                            "size": log_file_path.stat().st_size,
                            "modified_time": str(log_file_path.stat().st_mtime)
                        }
                    except Exception as e:
                        # 파일 읽기 실패 시에도 에러 정보 포함
                        manifest_data["execution_files"]["logs"][log_file_path.name] = {
                            "content": f"Error reading file: {str(e)}",
                            "size": log_file_path.stat().st_size if log_file_path.exists() else 0,
                            "modified_time": str(log_file_path.stat().st_mtime) if log_file_path.exists() else None,
                            "error": str(e)
                        }

        # Snakefile 내용 포함
        snakefile_path = os.path.join(task_folder_path, "Snakefile")
        if os.path.exists(snakefile_path):
            try:
                with open(snakefile_path, 'r', encoding='utf-8') as f:
                    manifest_data["execution_files"]["snakefile"] = {
                        "content": f.read(),
                        "path": snakefile_path
                    }
            except Exception as e:
                manifest_data["execution_files"]["snakefile"] = {
                    "content": f"Error reading Snakefile: {str(e)}",
                    "path": snakefile_path,
                    "error": str(e)
                }

        # Plugin metadata.json 파일 내용 포함
        plugin_metadata_path = os.path.join(task.plugin.plugin_path, "metadata.json")
        if os.path.exists(plugin_metadata_path):
            try:
                with open(plugin_metadata_path, 'r', encoding='utf-8') as f:
                    import json as json_lib
                    manifest_data["execution_files"]["plugin_metadata"] = {
                        "content": json_lib.load(f),
                        "path": plugin_metadata_path
                    }
            except Exception as e:
                manifest_data["execution_files"]["plugin_metadata"] = {
                    "content": f"Error reading plugin metadata.json: {str(e)}",
                    "path": plugin_metadata_path,
                    "error": str(e)
                }

        # meta.yml 파일 내용 포함 (있는 경우)
        # 아카이브된 로그를 사용하는 경우 아카이브된 meta.yml도 확인
        archived_meta_yml_path = os.path.join(task_folder_path, "executions", task.task_id, "meta.yml")
        if os.path.exists(archived_meta_yml_path):
            meta_yml_path = archived_meta_yml_path
        else:
            meta_yml_path = os.path.join(task_folder_path, "meta.yml")

        if os.path.exists(meta_yml_path):
            try:
                with open(meta_yml_path, 'r', encoding='utf-8') as f:
                    manifest_data["execution_files"]["meta_yml"] = {
                        "content": f.read(),
                        "path": meta_yml_path
                    }
            except Exception as e:
                manifest_data["execution_files"]["meta_yml"] = {
                    "content": f"Error reading meta.yml: {str(e)}",
                    "path": meta_yml_path,
                    "error": str(e)
                }

        # Results 디렉토리의 파일들 메타데이터 수집
        results_folder_path = os.path.join(task_folder_path, "results")
        if os.path.exists(results_folder_path) and os.path.isdir(results_folder_path):
            results_path = Path(results_folder_path)

            # 파일 타입 분류 함수
            def get_file_type(extension):
                file_type_mapping = {
                    '.sif': 'network',
                    '.txt': 'text',
                    '.csv': 'tabular',
                    '.tsv': 'tabular',
                    '.json': 'visualization',
                    '.npy': 'numpy_array',
                    '.h5ad': 'anndata',
                    '.png': 'image',
                    '.jpg': 'image',
                    '.jpeg': 'image',
                    '.pdf': 'document'
                }
                return file_type_mapping.get(extension.lower(), 'unknown')

            for result_file_path in results_path.glob("*"):
                if result_file_path.is_file():
                    try:
                        file_stat = result_file_path.stat()
                        file_extension = result_file_path.suffix

                        manifest_data["execution_files"]["results"][result_file_path.name] = {
                            "size": file_stat.st_size,
                            "modified_time": str(file_stat.st_mtime),
                            "extension": file_extension,
                            "file_type": get_file_type(file_extension),
                            "path": str(result_file_path)
                        }
                    except Exception as e:
                        # 파일 정보 읽기 실패 시에도 기본 정보는 포함
                        manifest_data["execution_files"]["results"][result_file_path.name] = {
                            "size": 0,
                            "modified_time": None,
                            "extension": result_file_path.suffix,
                            "file_type": "error",
                            "error": str(e)
                        }

            # 결과 파일이 없는 경우 정보 표시
            if not manifest_data["execution_files"]["results"]:
                logger.warning(f"No result files found in {results_folder_path}")
        else:
            logger.info(f"Results directory not found or not accessible: {results_folder_path}")

        # JSON으로 변환
        json_content = json.dumps(manifest_data, indent=2, ensure_ascii=False)

        # 파일명 생성 (타임스탬프 포함)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"execution_manifest_{task.plugin_name}_{task_id[:8]}_{timestamp}.json"

        # StreamingResponse로 다운로드 제공
        return StreamingResponse(
            io.StringIO(json_content),
            media_type="application/json",
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating execution manifest: {str(e)}")
