from fastapi import APIRouter, Depends, File, UploadFile, Form
from sqlalchemy.orm import Session
from typing import List
import logging

from app.auth import deps as dep
from app.plugin.schemas import PluginData, PluginCreate, PluginUpdate, PluginAssociate, BuildDockerRequest
from app import models
from app.plugin import service

# 로거 설정
logger = logging.getLogger(__name__)

# --- Test-patch compatibility re-exports ---------------------------------
# The business logic now lives in ``app.plugin.service``; the router only wires
# HTTP concerns to it. Existing tests (test_plugin_api.py) patch these symbols
# on the router namespace (e.g. ``app.plugin.router.AsyncResult``). Keeping them
# importable here preserves those patch targets without changing behavior.
# noqa: F401 imports below are intentional re-exports for test compatibility.
import os  # noqa: F401
import docker  # noqa: F401
from fastapi.responses import FileResponse  # noqa: F401
from celery.result import AsyncResult  # noqa: F401
from celery import current_app as celery_app  # noqa: F401
from app.plugin import utils as plugin_utils  # noqa: F401
from app.plugin import crud as crud_plugin  # noqa: F401
from app.plugin.utils import get_plugin_path, is_plugin_editable, ensure_local_plugins_dir  # noqa: F401
from app.plugin.cache import invalidate_all_plugin_cache  # noqa: F401
from app.worker.tasks import build_plugin_task  # noqa: F401

router = APIRouter()


@router.post("/validation")
def validate_plugin(
    *,
    plugin_data: PluginData,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.validate_plugin(plugin_data=plugin_data, username=current_user.username)


@router.post("/upload")
async def upload_plugin(
    *,
    db: Session = Depends(dep.get_db),
    plugin_data: PluginCreate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.upload_plugin(db=db, plugin_data=plugin_data)


@router.post("/upload_scripts")
async def upload_scripts(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return await service.upload_scripts(plugin_name=plugin_name, files=files)


@router.post("/upload_package")
async def upload_package(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
    ):
    return await service.upload_package(plugin_name=plugin_name, files=files)


@router.post("/build_docker/{plugin_name}")
async def build_plugin_docker(
    *,
    plugin_name: str,
    request: BuildDockerRequest,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    플러그인의 Dockerfile 생성 및 Docker 이미지 빌드
    스크립트와 패키지 파일들이 모두 업로드된 후에 실행되어야 함
    """
    return await service.build_plugin_docker(
        plugin_name=plugin_name,
        user_id=current_user.id,
        use_gpu=request.use_gpu,
    )


@router.get("/reference_folders/{plugin_name}")
def get_reference_folders(
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_reference_folders(plugin_name=plugin_name)


@router.get("/package/{plugin_name}")
def get_package_files(
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_package_files(plugin_name=plugin_name)


@router.get("/file/{plugin_name}/{file_name}")
def get_plugin_file(
    plugin_name: str,
    file_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_plugin_file(plugin_name=plugin_name, file_name=file_name)


@router.get("/list")
def list_plugins(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.list_plugins(db=db, user=current_user)


# plugin_name을 받아서 해당 플러그인의 정보를 반환
@router.get("/info/{plugin_name}")
def get_plugin_info(
    *,
    plugin_name: str,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_plugin_info(db=db, plugin_name=plugin_name)


@router.post("/associate")
def associate_plugin(
    *,
    db: Session = Depends(dep.get_db),
    pluginInfo: PluginAssociate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.associate_plugin(db=db, user=current_user, plugin_id=pluginInfo.plugin_id)


@router.post("/dissociate")
def dissociate_plugin(
    *,
    db: Session = Depends(dep.get_db),
    pluginInfo: PluginAssociate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.dissociate_plugin(db=db, user=current_user, plugin_id=pluginInfo.plugin_id)


@router.get("/template/{plugin_id}")
def get_plugin_template(
    *,
    db: Session = Depends(dep.get_db),
    plugin_id: int,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return service.get_plugin_template(db=db, plugin_id=plugin_id)


@router.post("/build/{plugin_name}")
async def build_plugin(
    *,
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    플러그인 Docker 이미지를 비동기적으로 빌드합니다.
    Celery task를 사용하여 백그라운드에서 처리됩니다.
    """
    return await service.build_plugin(plugin_name=plugin_name, user_id=current_user.id)


@router.get("/check_image/{plugin_name}")
async def check_plugin_image(
    *,
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    return await service.check_plugin_image(plugin_name=plugin_name)


@router.post("/update/{plugin_name}")
async def update_plugin_complete(
    *,
    plugin_name: str,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    플러그인의 메타데이터를 물리적 파일에서 DB로 동기화하는 엔드포인트
    스크립트나 패키지 파일 업데이트 후 DB 정보를 최신화할 때 사용
    """
    return await service.update_plugin_complete(db=db, plugin_name=plugin_name)


@router.post("/upload_text_dependencies")
async def upload_text_dependencies(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
    db: Session = Depends(dep.get_db),
):
    """
    텍스트 의존성 파일들(requirements.txt, environment.yml, renv.lock)을 업로드하고 DB도 업데이트
    """
    return await service.upload_text_dependencies(db=db, plugin_name=plugin_name, files=files)


@router.get("/build/status/{task_id}")
async def get_build_status(
    *,
    task_id: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    플러그인 빌드 태스크의 상태를 조회합니다.
    """
    return await service.get_build_status(task_id=task_id)


@router.get("/build/tasks")
async def get_build_tasks(
    *,
    plugin_name: str = None,
    current_user: models.User = Depends(dep.get_current_active_user),
    db: Session = Depends(dep.get_db),
):
    """
    플러그인 빌드 태스크 목록을 조회합니다.
    plugin_name이 제공되면 해당 플러그인의 태스크만 필터링합니다.
    """
    return await service.get_build_tasks(db=db, user=current_user, plugin_name=plugin_name)


@router.post("/build/cancel/{task_id}")
async def cancel_build_task(
    *,
    task_id: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    실행 중인 플러그인 빌드 태스크를 취소합니다.
    """
    return await service.cancel_build_task(task_id=task_id)


@router.get("/build/logs/{plugin_name}")
async def get_build_logs(
    *,
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    """
    플러그인 빌드 로그를 조회합니다.
    """
    return await service.get_build_logs(plugin_name=plugin_name)


# Version management endpoints for official plugins
@router.get("/versions/{plugin_name}")
async def get_plugin_versions(
    *,
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
    db: Session = Depends(dep.get_db),
):
    """
    Get available versions for an official plugin from GitHub Container Registry.
    """
    return await service.get_plugin_versions(db=db, plugin_name=plugin_name)


@router.post("/versions/{plugin_name}/update")
async def update_plugin_version(
    *,
    plugin_name: str,
    version: str = Form(...),
    current_user: models.User = Depends(dep.get_current_active_user),
    db: Session = Depends(dep.get_db),
):
    """
    Update the version of an official plugin.
    """
    return await service.update_plugin_version(db=db, plugin_name=plugin_name, version=version)


@router.get("/version/{plugin_name}")
async def get_plugin_current_version(
    *,
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
    db: Session = Depends(dep.get_db),
):
    """
    Get current version information for a plugin.
    """
    return await service.get_plugin_current_version(db=db, plugin_name=plugin_name)
