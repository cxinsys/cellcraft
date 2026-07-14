"""
Plugin domain service layer.

Extracted from ``plugin/router.py`` in PR-5 (Phase 3a). This module holds the
business logic that previously lived inline in the endpoints:
filesystem I/O, metadata/Snakefile generation, DB upserts + rollback, Docker
build orchestration (delegated to ``plugin/builder.py``), and version
management. The router now only parses requests, resolves dependencies, calls a
service function, and returns its result.

Transitional exception policy (PR-5): the global domain-exception layer arrives
in PR-8. Until then, service functions may ``raise HTTPException`` directly
(copy-move fidelity first, exception translation later). Response shapes — JSON
keys, status codes, and error message strings — are preserved exactly.
"""
import os
import time
import json
import shutil
import pathlib
import uuid
import logging
from datetime import datetime
from typing import List, Optional

from fastapi import HTTPException, UploadFile
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from celery.result import AsyncResult
from celery import current_app as celery_app

from app import models
from app.plugin import crud as crud_plugin
from app.plugin import utils as plugin_utils
from app.plugin import builder
from app.plugin.utils import get_plugin_path, is_plugin_editable, ensure_local_plugins_dir
from app.plugin.cache import invalidate_all_plugin_cache
from app.plugin.schemas import PluginData, PluginCreate

logger = logging.getLogger(__name__)


def parse_reference_folders(folders):
    """
    재귀적으로 referenceFolders 데이터를 dict 형식으로 변환.

    Parameters:
        folders (list): ReferenceFolders 리스트.

    Returns:
        dict: 변환된 딕셔너리.
    """
    result = {}

    for folder in folders:
        folder_dict = {
            file.fileName: file.file for file in folder.files  # 파일 저장
        }

        # 하위 폴더 리스트를 새롭게 생성하여 중복 방지
        sub_folders = []
        for subFolder in folder.subFolders:
            sub_folders.append(parse_reference_folders([subFolder]))  # 새로운 리스트에 저장

        folder_dict["subFolders"] = sub_folders  # subFolders 속성 추가
        result[folder.folderName] = folder_dict  # 폴더 이름을 키로 추가

    return result


def validate_plugin(*, plugin_data: PluginData, username: str) -> dict:
    """Validate plugin metadata and produce a PluginCreate draft (as a local plugin)."""
    try:
        # Convert PluginInfo.dependencyFiles to a dictionary
        dependencies_dict = {dep.fileName: dep.file for dep in plugin_data.plugin.dependencyFiles}
        reference_folders = parse_reference_folders(plugin_data.plugin.referenceFolders)
        rules_dict = {
            index: {
                "name": rule.name,
                "input": rule.input,
                "output": rule.output,
                "script": rule.script,
                "parameters": rule.parameters,
                "nodeId": rule.nodeId,
                "isVisualization": rule.isVisualization
            }
            for index, rule in enumerate(plugin_data.rules)
        }

        # 기본 검증
        if not plugin_data.plugin.name or not plugin_data.plugin.description:
            raise HTTPException(
                status_code=400,
                detail="Plugin name and description are required"
            )

        # 규칙 검증
        if not plugin_data.rules:
            raise HTTPException(
                status_code=400,
                detail="At least one rule is required"
            )

        # 스크립트 파일 검증
        scripts = [rule.script for rule in plugin_data.rules if rule.script]
        if not scripts:
            raise HTTPException(
                status_code=400,
                detail="At least one script file is required"
            )

        # Convert PluginData to PluginCreate - always create as local plugin
        plugin_upload_data = PluginCreate(
            name=plugin_data.plugin.name,
            description=plugin_data.plugin.description,
            author=username,
            plugin_path=f"./plugin/local/{plugin_data.plugin.name}/",
            plugin_type=plugin_data.plugin.pluginType if hasattr(plugin_data.plugin, 'pluginType') else None,
            dependencies=dependencies_dict if dependencies_dict else None,
            reference_folders=reference_folders if reference_folders else None,
            drawflow=plugin_data.drawflow,
            rules=rules_dict,
            use_gpu=plugin_data.plugin.useGpu if hasattr(plugin_data.plugin, 'useGpu') else False,
            source="local",
            is_editable=True
        )

        return {
            "message": "Plugin data validated",
            "plugin": plugin_upload_data,
            "scripts": scripts
        }
    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Validation error: {str(e)}")


# ---------------------------------------------------------------------------
# upload_plugin: decomposed into step functions (file save / metadata+snakefile
# / DB upsert) with an explicit rollback path.
# ---------------------------------------------------------------------------

def _backup_existing_plugin(*, plugin_data: PluginCreate, plugin_folder: str,
                            db_existing_plugin) -> Optional[str]:
    """Back up an existing local plugin folder before recreating it.

    Returns the backup folder path (or ``None`` when no backup was made). On
    backup failure the partial backup is cleaned up and an HTTPException(500) is
    raised — matching the original inline behavior.
    """
    if db_existing_plugin and os.path.exists(plugin_folder):
        # 상위 디렉토리에 백업 폴더 생성
        backup_base = "./plugin/backups"
        if not os.path.exists(backup_base):
            os.makedirs(backup_base)

        backup_folder = os.path.join(backup_base, f"{plugin_data.name}_backup_{int(time.time())}")
        try:
            shutil.copytree(plugin_folder, backup_folder)
            shutil.rmtree(plugin_folder)
        except Exception as e:
            if os.path.exists(backup_folder):
                shutil.rmtree(backup_folder)
            raise HTTPException(
                status_code=500,
                detail=f"기존 플러그인 백업 실패: {str(e)}"
            )
        return backup_folder
    return None


def _save_plugin_files(*, plugin_data: PluginCreate, plugin_folder: str,
                       dependency_folder: str, script_folder: str,
                       backup_folder: Optional[str]) -> None:
    """Create plugin/dependency/scripts folders and reference folders.

    On failure, restores from ``backup_folder`` if present and raises
    HTTPException(500) — matching the original inline behavior.
    """
    try:
        # 3. 새 폴더 생성
        print(f"Creating plugin folder: {plugin_folder}")
        plugin_utils.create_plugin_folder(plugin_folder)

        print(f"Creating dependency folder: {dependency_folder}")
        plugin_utils.create_dependency_folder(dependency_folder, plugin_data.dependencies)

        # scripts 폴더는 항상 생성 (Docker 빌드를 위해 필요)
        print(f"Ensuring scripts folder exists: {script_folder}")
        if not os.path.exists(script_folder):
            os.makedirs(script_folder)
            print(f"Created scripts folder: {script_folder}")

        if plugin_data.reference_folders:
            print(f"Creating reference folders in: {script_folder}")
            plugin_utils.create_reference_folder(script_folder, plugin_data.reference_folders)
        else:
            print("No reference folders provided")
    except Exception as e:
        # 폴더 생성 실패 시 백업 복원
        if backup_folder and os.path.exists(backup_folder):
            if os.path.exists(plugin_folder):
                shutil.rmtree(plugin_folder)
            shutil.copytree(backup_folder, plugin_folder)
            shutil.rmtree(backup_folder)
        raise HTTPException(
            status_code=500,
            detail=f"플러그인 폴더 생성 실패: {str(e)}"
        )


def _rollback_upload(*, plugin_folder: str, backup_folder: Optional[str]) -> None:
    """Restore the previous plugin state after a failed metadata/DB step.

    Removes the freshly-created plugin folder and, if a backup exists, restores
    it in place — matching the inner ``except Exception`` rollback of the
    original ``upload_plugin``.
    """
    if os.path.exists(plugin_folder):
        shutil.rmtree(plugin_folder)
    if backup_folder and os.path.exists(backup_folder):
        shutil.copytree(backup_folder, plugin_folder)
        shutil.rmtree(backup_folder)


def _upsert_plugin(*, db: Session, plugin_data: PluginCreate, plugin_folder: str,
                   db_existing_plugin, backup_folder: Optional[str]) -> dict:
    """Write metadata + Snakefile, upsert the DB row, commit, and clean up backup.

    On failure rolls back the DB transaction and the filesystem via
    ``_rollback_upload`` and raises HTTPException(500). Returns the success
    response dict on success.
    """
    try:
        # 4. 메타데이터 생성
        metadata = {
            "name": plugin_data.name,
            "author": plugin_data.author,
            "description": plugin_data.description,
            "drawflow": plugin_data.drawflow,
            "rules": plugin_data.rules
        }
        plugin_utils.create_metadata_file(plugin_folder, metadata)

        # 5. Snakefile 생성
        # For plugin creation, we don't have specific workflow/visualization IDs,
        # so we generate a generic Snakefile that will be customized during execution
        plugin_utils.generate_snakemake_code(
            plugin_data.rules,
            plugin_folder,
            plugin_type=plugin_data.plugin_type
        )

        # 6. 데이터베이스 업데이트
        if db_existing_plugin:
            db_plugin = crud_plugin.update_plugin(db=db, plugin=plugin_data, plugin_id=db_existing_plugin.id)
        else:
            db_plugin = crud_plugin.create_plugin(db=db, plugin=plugin_data)

        # 7. 트랜잭션 커밋
        db.commit()

        # 8. 성공 시 백업 삭제
        if backup_folder and os.path.exists(backup_folder):
            shutil.rmtree(backup_folder)

        invalidate_all_plugin_cache()
        return {
            "message": "플러그인 메타데이터 업로드 성공",
            "plugin": db_plugin,
            "success": True,
            "redirect_to": "plugin_list"  # 프론트엔드에서 플러그인 페이지로 리다이렉트하도록 지시
        }

    except HTTPException as he:
        # HTTP 예외는 그대로 전달
        raise he
    except Exception as e:
        # 실패 시 롤백
        db.rollback()
        _rollback_upload(plugin_folder=plugin_folder, backup_folder=backup_folder)

        raise HTTPException(status_code=500, detail=f"플러그인 메타데이터 업데이트 실패: {str(e)}")


def upload_plugin(*, db: Session, plugin_data: PluginCreate) -> dict:
    """Upload (create or update) a local plugin: files + metadata + DB, with rollback."""
    plugin_folder = None
    backup_folder = None
    try:
        # 1. Ensure local plugins directory exists and create plugin folder
        ensure_local_plugins_dir()

        # Check if plugin already exists and is editable
        try:
            existing_path, existing_source = get_plugin_path(plugin_data.name)
            if existing_source == "official":
                raise HTTPException(
                    status_code=403,
                    detail=f"Cannot modify official plugin '{plugin_data.name}'. Official plugins are read-only."
                )
        except HTTPException:
            # Plugin doesn't exist yet, which is fine for upload
            pass

        plugin_folder = f"./plugin/local/{plugin_data.name}/"
        dependency_folder = os.path.join(plugin_folder, "dependency")
        script_folder = os.path.join(plugin_folder, "scripts")

        # 2. 기존 플러그인 확인 및 백업 (source를 고려하여 검색)
        db_existing_plugin = db.query(models.Plugin).filter(
            models.Plugin.name == plugin_data.name,
            models.Plugin.source == "local"  # 로컬 플러그인만 업데이트 가능
        ).first()

        # Official 플러그인과 이름이 같은지 확인
        official_plugin = db.query(models.Plugin).filter(
            models.Plugin.name == plugin_data.name,
            models.Plugin.source == "official"
        ).first()

        if official_plugin:
            # Official 플러그인과 같은 이름의 로컬 플러그인은 허용하지만 경고
            print(f"Warning: Creating local plugin with same name as official plugin: {plugin_data.name}")

        backup_folder = _backup_existing_plugin(
            plugin_data=plugin_data,
            plugin_folder=plugin_folder,
            db_existing_plugin=db_existing_plugin,
        )

        _save_plugin_files(
            plugin_data=plugin_data,
            plugin_folder=plugin_folder,
            dependency_folder=dependency_folder,
            script_folder=script_folder,
            backup_folder=backup_folder,
        )

        return _upsert_plugin(
            db=db,
            plugin_data=plugin_data,
            plugin_folder=plugin_folder,
            db_existing_plugin=db_existing_plugin,
            backup_folder=backup_folder,
        )

    except HTTPException as he:
        raise he
    except Exception as e:
        if plugin_folder and os.path.exists(plugin_folder):
            shutil.rmtree(plugin_folder)
        if backup_folder and os.path.exists(backup_folder):
            shutil.rmtree(backup_folder)

        raise HTTPException(status_code=500, detail=f"플러그인 업로드 중 예기치 않은 오류: {str(e)}")


async def upload_scripts(*, plugin_name: str, files: List[UploadFile]) -> dict:
    """Stage, promote, and (on failure) roll back a plugin scripts upload."""
    # Check if plugin is editable
    if not is_plugin_editable(plugin_name):
        raise HTTPException(
            status_code=403,
            detail=f"Cannot modify plugin '{plugin_name}'. Official plugins are read-only."
        )

    # Get the actual plugin path (should be local)
    plugin_path_str, source = get_plugin_path(plugin_name, "local")
    plugin_path = pathlib.Path(plugin_path_str)

    # 디렉토리 이름 정의
    final_scripts_dirname = "scripts"  # 최종 스크립트가 위치할 디렉토리 이름
    previous_scripts_dirname = "scripts_previous"
    staging_scripts_basename = "scripts_staging"

    # 전체 경로 정의
    scripts_dir = plugin_path / final_scripts_dirname
    scripts_previous_dir = plugin_path / previous_scripts_dirname
    # 각 업로드 시도마다 고유한 스테이징 디렉토리 생성
    scripts_staging_dir = plugin_path / f"{staging_scripts_basename}_{uuid.uuid4()}"

    try:
        logger.info(f"Starting script upload for plugin: {plugin_name}. Target directory: {scripts_dir}")

        # 1. 스테이징 디렉토리 생성
        scripts_staging_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created staging directory: {scripts_staging_dir}")

        # 2. 스크립트 파일을 스테이징 디렉토리에 저장
        for file_obj in files:
            file_path = scripts_staging_dir / file_obj.filename
            try:
                # UploadFile.read()는 async 메소드이므로 await 사용
                content = await file_obj.read()
                with open(file_path, "wb") as f:
                    f.write(content)
                logger.info(f"Successfully saved file {file_obj.filename} to {scripts_staging_dir}")
            except Exception as e:
                logger.error(f"Failed to save file {file_obj.filename} to staging directory {scripts_staging_dir}: {str(e)}")
                if scripts_staging_dir.exists():
                    shutil.rmtree(scripts_staging_dir)
                    logger.info(f"Cleaned up staging directory {scripts_staging_dir} due to file save error.")
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to save script file {file_obj.filename}: {str(e)}"
                )
        logger.info(f"All files successfully saved to staging directory: {scripts_staging_dir}")

        # 3. 기존 스크립트 백업 (scripts_dir -> scripts_previous_dir)
        if scripts_dir.exists():
            logger.info(f"Current scripts directory {scripts_dir} exists. Proceeding with backup.")
            if scripts_previous_dir.exists():
                logger.info(f"Removing old previous scripts directory: {scripts_previous_dir}")
                shutil.rmtree(scripts_previous_dir)

            # 전체 scripts_dir를 scripts_previous_dir로 이동 (rename)
            logger.info(f"Moving entire {scripts_dir} to {scripts_previous_dir}")
            shutil.move(str(scripts_dir), str(scripts_previous_dir))
            logger.info(f"Backup of {scripts_dir} to {scripts_previous_dir} completed.")
        else:
            logger.info(f"No current scripts directory found at {scripts_dir}. Skipping backup step.")

        # --- 핵심 트랜잭션: 스테이징 디렉토리를 최종 디렉토리로 승격 ---
        try:
            # 4. scripts_dir이 존재하지 않는 상태에서 스테이징 디렉토리를 최종 스크립트 디렉토리로 이동
            logger.info(f"Moving staging scripts {scripts_staging_dir} to final target {scripts_dir}")

            # 이제 scripts_dir이 존재하지 않으므로 rename 동작이 수행됨
            shutil.move(str(scripts_staging_dir), str(scripts_dir))

            # 5. (선택사항) 기존 scripts_previous에서 필요한 하위 디렉토리(예: reference) 복원
            if scripts_previous_dir.exists():
                for item in scripts_previous_dir.iterdir():
                    if item.is_dir():  # 디렉토리인 경우
                        target_dir = scripts_dir / item.name
                        if not target_dir.exists():  # 새 scripts에 없는 디렉토리만 복원
                            logger.info(f"Restoring directory {item.name} from previous scripts")
                            shutil.copytree(str(item), str(target_dir))

            logger.info(f"Successfully moved {scripts_staging_dir} to {scripts_dir}. Script update complete.")

            invalidate_all_plugin_cache()
            return {
                "message": "Scripts uploaded successfully",
                "scripts_path": str(scripts_dir),
                "scripts": [file.filename for file in files],
                "success": True
            }

        except Exception as e_promote:
            logger.error(f"Failed to move staging directory {scripts_staging_dir} to {scripts_dir}: {str(e_promote)}")

            # 롤백 시도
            logger.info("Attempting rollback...")
            try:
                # 부분적으로 생성되었을 수 있는 새 scripts_dir 삭제
                if scripts_dir.exists():
                    logger.info(f"Removing potentially incomplete/failed {scripts_dir}")
                    shutil.rmtree(scripts_dir)

                # scripts_previous_dir가 존재하면 (백업이 있었다면) 원래 위치로 복원
                if scripts_previous_dir.exists():
                    logger.info(f"Attempting to restore {scripts_previous_dir} to {scripts_dir}")
                    shutil.move(str(scripts_previous_dir), str(scripts_dir))
                    logger.info(f"Rollback successful: {scripts_previous_dir} restored to {scripts_dir}.")
                else:
                    logger.warning(f"Rollback attempted, but no previous scripts directory ({scripts_previous_dir}) found to restore.")

            except Exception as e_rollback:
                logger.critical(f"Rollback failed: {str(e_rollback)}. Original promotion error: {str(e_promote)}")
                # 롤백 실패 시, 관리자 알림 등의 로직 추가 고려
                if scripts_staging_dir.exists():
                    shutil.rmtree(scripts_staging_dir)
                    logger.info(f"Cleaned up staging directory {scripts_staging_dir} after failed promotion and failed rollback.")
                raise HTTPException(
                    status_code=500,
                    detail={
                        "message": "Failed to update scripts, and subsequent rollback also failed.",
                        "promotion_error": str(e_promote),
                        "rollback_error": str(e_rollback)
                    }
                )

            # 롤백 시도 후 (성공 여부와 관계없이) 스테이징 디렉토리는 실패한 버전이므로 정리
            if scripts_staging_dir.exists():
                shutil.rmtree(scripts_staging_dir)
                logger.info(f"Cleaned up staging directory {scripts_staging_dir} after failed promotion.")

            raise HTTPException(
                status_code=500,
                detail=f"Failed to update scripts. Rollback attempted. Original error: {str(e_promote)}"
            )
    # --- 핵심 트랜잭션 종료 ---

    except HTTPException as he:
        # 이미 로깅된 HTTP 예외가 아니라면 로깅
        logger.warning(f"HTTPException caught during script upload: {he.detail}")
        # 스테이징 디렉토리가 남아있다면 정리
        if scripts_staging_dir.exists():
            shutil.rmtree(scripts_staging_dir)
            logger.info(f"Cleaned up staging directory {scripts_staging_dir} due to HTTPException.")
        raise he

    except Exception as e_unexpected:
        logger.critical(f"Unexpected error during script upload for plugin {plugin_name}: {str(e_unexpected)}", exc_info=True)
        # 예상치 못한 오류 발생 시, 최소한 스테이징 디렉토리는 정리
        if scripts_staging_dir.exists():
            try:
                shutil.rmtree(scripts_staging_dir)
                logger.info(f"Cleaned up staging directory {scripts_staging_dir} due to unexpected error.")
            except Exception as e_cleanup_staging:
                logger.error(f"Failed to cleanup staging directory {scripts_staging_dir} during unexpected error handling: {str(e_cleanup_staging)}")

        raise HTTPException(
            status_code=500,
            detail=f"An unexpected error occurred during script upload: {str(e_unexpected)}"
        )


async def upload_package(*, plugin_name: str, files: List[UploadFile]) -> dict:
    """Upload package (.whl/.tar.gz) files into a plugin's dependency folder."""
    try:
        # Check if plugin is editable
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Cannot modify plugin '{plugin_name}'. Official plugins are read-only."
            )

        # Get the actual plugin path (should be local)
        plugin_path_str, source = get_plugin_path(plugin_name, "local")
        dependency_folder = os.path.join(plugin_path_str, "dependency")
        temp_folder = os.path.join(plugin_path_str, "dependency_temp")
        backup_folder = None

        # 1. dependency 폴더가 존재하지 않으면 생성
        if not os.path.exists(dependency_folder):
            os.makedirs(dependency_folder)

        # 2. 임시 폴더 생성
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

        # 3. 패키지 파일을 임시 폴더에 저장
        uploaded_file_names = []
        for file in files:
            if not file.filename.endswith(('.whl', '.tar.gz')):
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid package file format: {file.filename}"
                )

            uploaded_file_names.append(file.filename)
            file_path = os.path.join(temp_folder, file.filename)
            try:
                with open(file_path, "wb") as f:
                    content = await file.read()
                    f.write(content)
            except Exception as e:
                # 실패 시 임시 폴더 삭제
                if os.path.exists(temp_folder):
                    shutil.rmtree(temp_folder)
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to save package file {file.filename}: {str(e)}"
                )

        # 4. 기존 dependency 폴더 백업
        if os.path.exists(dependency_folder):
            backup_folder = f"{dependency_folder}_backup_{int(time.time())}"
            shutil.copytree(dependency_folder, backup_folder)

        try:
            # 5. 기존 dependency 폴더에서 동일한 이름의 패키지 파일만 제거 (선택적 교체)
            if os.path.exists(dependency_folder):
                for uploaded_filename in uploaded_file_names:
                    existing_file_path = os.path.join(dependency_folder, uploaded_filename)
                    if os.path.isfile(existing_file_path):
                        os.remove(existing_file_path)
                        print(f"Removed existing package file for replacement: {uploaded_filename}")

            # 6. 새로운 패키지 파일들을 dependency 폴더로 복사
            for file in files:
                src_path = os.path.join(temp_folder, file.filename)
                dst_path = os.path.join(dependency_folder, file.filename)
                shutil.copy2(src_path, dst_path)
                print(f"Copied new package file: {file.filename}")

            # 7. 임시 폴더 삭제
            if os.path.exists(temp_folder):
                shutil.rmtree(temp_folder)

            # 8. 성공 시 백업 폴더 삭제
            if backup_folder and os.path.exists(backup_folder):
                shutil.rmtree(backup_folder)

            # 최종 파일 목록 출력 (디버깅용)
            final_package_files = [f for f in os.listdir(dependency_folder) if f.endswith(('.whl', '.tar.gz'))]
            print(f"Final package files in dependency folder: {final_package_files}")

            invalidate_all_plugin_cache()
            return {
                "message": "Package uploaded successfully",
                "packages": [file.filename for file in files],
                "success": True
            }

        except Exception as e:
            # 실패 시 복구
            if backup_folder and os.path.exists(backup_folder):
                # 기존 dependency 폴더 삭제 후 백업에서 복원
                if os.path.exists(dependency_folder):
                    shutil.rmtree(dependency_folder)
                shutil.copytree(backup_folder, dependency_folder)
                shutil.rmtree(backup_folder)

            # 임시 폴더 정리
            if os.path.exists(temp_folder):
                shutil.rmtree(temp_folder)

            raise HTTPException(
                status_code=500,
                detail=f"Failed to update package files: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Unexpected error during package upload: {str(e)}"
        )


async def build_plugin_docker(*, plugin_name: str, user_id: int, use_gpu: bool) -> dict:
    """Generate the plugin Dockerfile and start an async Docker image build."""
    try:
        # Check if plugin is editable for build operations
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Cannot build official plugin '{plugin_name}'. Only local plugins can be built."
            )

        # Get the actual plugin path (should be local)
        plugin_folder, source = get_plugin_path(plugin_name, "local")
        script_folder = os.path.join(plugin_folder, "scripts")

        # 플러그인 폴더가 존재하는지 확인
        if not os.path.exists(plugin_folder):
            raise HTTPException(
                status_code=404,
                detail=f"Plugin folder not found: {plugin_name}"
            )

        # scripts 폴더 준비 및 Dockerfile 생성 (builder에 위임)
        builder.prepare_build_context(
            plugin_folder=plugin_folder,
            script_folder=script_folder,
            use_gpu=use_gpu,
        )

        # 비동기 빌드 태스크 시작
        task = builder.dispatch_build_task(plugin_name=plugin_name, user_id=user_id)

        invalidate_all_plugin_cache()
        return {
            "message": f"플러그인 Docker 이미지 빌드 시작: {plugin_name}",
            "task_id": task.id,
            "plugin_name": plugin_name,
            "status": "PENDING"
        }

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to build plugin Docker image: {str(e)}"
        )


def get_reference_folders(*, plugin_name: str) -> dict:
    """List a plugin's reference folder structures under scripts/."""
    try:
        # Get plugin path (can be either official or local for reading)
        plugin_path, source = get_plugin_path(plugin_name)
        folder_path = os.path.join(plugin_path, "scripts")

        # 하위 폴더 리스트 가져오기
        folder_names = plugin_utils.get_reference_folders_list(folder_path)

        # 각 폴더의 구조 가져오기
        reference_folders = []
        for folder_name in folder_names:
            sub_folder_path = os.path.join(folder_path, folder_name)
            folder_structure = plugin_utils.get_reference_folder(sub_folder_path)
            reference_folders.append(folder_structure)

        return {"reference_folders": reference_folders}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def get_package_files(*, plugin_name: str) -> dict:
    """List package files (.whl/.tar.gz) in a plugin's dependency folder."""
    try:
        # Get plugin path (can be either official or local for reading)
        plugin_path, source = get_plugin_path(plugin_name)
        folder_path = os.path.join(plugin_path, "dependency")

        # 파일 리스트 중에서 .whl, .tar.gz 파일만 가져와서 반환
        package_files = [file for file in os.listdir(folder_path) if file.endswith((".whl", ".tar.gz"))]
        return {"package_files": package_files}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def get_plugin_file(*, plugin_name: str, file_name: str) -> FileResponse:
    """Resolve and return a FileResponse for a file within a plugin's folders."""
    try:
        # Get plugin path (can be either official or local for reading)
        base_folder_path, source = get_plugin_path(plugin_name)

        # 폴더 안에 dependency 폴더, scripts 폴더
        dependency_folder_path = os.path.join(base_folder_path, "dependency")
        scripts_folder_path = os.path.join(base_folder_path, "scripts")

        # os.path.isfile(file_path)를 사용해서 각 폴더 내에 file_name에 해당하는 파일이 있는지 확인 후, 파일 경로 설정
        file_path = None
        if os.path.isfile(os.path.join(base_folder_path, file_name)):
            file_path = os.path.join(base_folder_path, file_name)
        elif os.path.isfile(os.path.join(dependency_folder_path, file_name)):
            file_path = os.path.join(dependency_folder_path, file_name)
        elif os.path.isfile(os.path.join(scripts_folder_path, file_name)):
            file_path = os.path.join(scripts_folder_path, file_name)
        else:
            # scripts 폴더의 하위 폴더를 재귀적으로 탐색하여 파일 찾기
            for root, _, files in os.walk(scripts_folder_path):
                if file_name in files:
                    file_path = os.path.join(root, file_name)
                    break

        # 파일 경로를 찾지 못했을 경우 예외 처리
        if not file_path:
            raise HTTPException(status_code=404, detail="File not found")

        return FileResponse(file_path, filename=file_name)

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def list_plugins(*, db: Session, user: models.User) -> dict:
    """Return the (cached) plugin list for a user, enriched with build/version info."""
    try:
        # Cache hit path
        from app.plugin.cache import get_cached_plugin_list, set_cached_plugin_list
        cached = get_cached_plugin_list(user.id)
        if cached is not None:
            return cached

        plugins = crud_plugin.get_plugins(db)
        plugin_list = []

        # Query user tasks once outside the loop (was N times inside)
        from app.task.crud import get_user_task
        try:
            user_tasks = get_user_task(db, user.id)
        except Exception as e:
            logger.warning(f"Error fetching user tasks: {e}")
            user_tasks = []

        # Create Docker client only if local plugins exist (avoids 21ms socket overhead)
        import docker
        docker_client = None
        has_local = any(p.source == "local" for p in plugins)
        if has_local:
            try:
                docker_client = docker.from_env()
            except Exception:
                docker_client = None

        for plugin in plugins:
            plugin_dict = plugin.__dict__
            plugin_dict['users'] = [user.__dict__ for user in plugin.users]

            # Convert rules dictionary to an array
            if 'rules' in plugin_dict:
                rules_dict = plugin_dict['rules']
                rules_array = [rules_dict[str(i)] for i in range(len(rules_dict))]
                plugin_dict['rules'] = rules_array

            # 빌드 상태 정보 추가 (only for local plugins)
            if plugin.source == "local" and docker_client:
                try:
                    image_tag = f"plugin-{plugin.name.lower()}"
                    docker_client.images.get(image_tag)
                    plugin_dict['docker_image_exists'] = True
                except docker.errors.ImageNotFound:
                    plugin_dict['docker_image_exists'] = False
                except Exception as e:
                    logger.warning(f"Error checking Docker image for {plugin.name}: {e}")
                    plugin_dict['docker_image_exists'] = False
            else:
                plugin_dict['docker_image_exists'] = False

            # 최근 빌드 태스크 정보 조회 (using pre-fetched user_tasks)
            try:
                plugin_build_tasks = [
                    task for task in user_tasks
                    if task.task_type == "plugin_build" and task.plugin_name == plugin.name
                ]

                if plugin_build_tasks:
                    latest_task = max(plugin_build_tasks, key=lambda x: x.start_time or datetime.min.replace(tzinfo=None))
                    task_result = AsyncResult(latest_task.task_id)

                    plugin_dict['latest_build'] = {
                        "task_id": latest_task.task_id,
                        "status": task_result.state,
                        "start_time": latest_task.start_time.isoformat() if latest_task.start_time else None,
                        "end_time": latest_task.end_time.isoformat() if latest_task.end_time else None
                    }
                else:
                    plugin_dict['latest_build'] = None
            except Exception as e:
                logger.warning(f"Error getting build tasks for {plugin.name}: {e}")
                plugin_dict['latest_build'] = None

            # Official 플러그인의 경우 데이터베이스 version 컬럼 직접 사용
            if plugin.source == "official":
                plugin_dict['current_version'] = plugin.version or "1.0"
                plugin_dict['available_versions'] = []

                logger.info(f"Plugin {plugin.name} - Using database version: {plugin_dict['current_version']}")
            else:
                plugin_dict['current_version'] = plugin.version or "local"
                plugin_dict['available_versions'] = []

            plugin_list.append(plugin_dict)

        if docker_client:
            docker_client.close()

        response = {"plugins": plugin_list}
        set_cached_plugin_list(user.id, response)
        return response
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def get_plugin_info(*, db: Session, plugin_name: str) -> dict:
    """Return a plugin's DB record by name (404 if missing)."""
    try:
        # Get the plugin by name
        plugin = crud_plugin.get_plugin_by_name(db, plugin_name)

        if plugin is None:
            raise HTTPException(status_code=404, detail="Plugin not found")

        return {"plugin": plugin}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def associate_plugin(*, db: Session, user: models.User, plugin_id: int) -> dict:
    """Associate a plugin with the current user."""
    try:
        # Associate the plugin with the current user
        crud_plugin.associate_user_plugin(db, user.id, plugin_id)

        # Get the plugin by ID
        plugin = crud_plugin.get_plugin_by_id(db, plugin_id)

        invalidate_all_plugin_cache()
        return {"message": "Plugin associated with user successfully", "plugin": plugin}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def dissociate_plugin(*, db: Session, user: models.User, plugin_id: int) -> dict:
    """Dissociate a plugin from the current user."""
    try:
        # Dissociate the plugin with the current user
        crud_plugin.dissociate_user_plugin(db, user.id, plugin_id)

        # Get the plugin by ID
        plugin = crud_plugin.get_plugin_by_id(db, plugin_id)

        invalidate_all_plugin_cache()
        return {"message": "Plugin dissociated from user successfully", "plugin": plugin}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def get_plugin_template(*, db: Session, plugin_id: int) -> dict:
    """Generate a drawflow template for a plugin (404 if missing)."""
    try:
        # Get the plugin by ID
        plugin = crud_plugin.get_plugin_by_id(db, plugin_id)

        if plugin is None:
            raise HTTPException(status_code=404, detail="Plugin not found")

        # Convert plugin.drawflow to dictionary if it's a JSON string
        if isinstance(plugin.drawflow, str):
            plugin.drawflow = json.loads(plugin.drawflow)

        # Generate the drawflow template
        drawflow = plugin_utils.generate_plugin_drawflow_template(plugin.drawflow, plugin.name)

        return {"drawflow": drawflow}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


async def build_plugin(*, plugin_name: str, user_id: int) -> dict:
    """Start an async Docker build for a plugin whose Dockerfile already exists."""
    try:
        # Check if plugin is editable for build operations
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Cannot build official plugin '{plugin_name}'. Only local plugins can be built."
            )

        # Get the actual plugin path (should be local)
        plugin_folder, source = get_plugin_path(plugin_name, "local")

        # 플러그인 폴더가 존재하는지 확인
        if not os.path.exists(plugin_folder):
            raise HTTPException(
                status_code=404,
                detail=f"Plugin folder not found: {plugin_name}"
            )

        # Dockerfile이 존재하는지 확인
        dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
        if not os.path.exists(dockerfile_path):
            raise HTTPException(
                status_code=404,
                detail=f"Dockerfile not found in plugin folder: {plugin_name}"
            )

        # 비동기 빌드 태스크 시작
        task = builder.dispatch_build_task(plugin_name=plugin_name, user_id=user_id)

        invalidate_all_plugin_cache()
        return {
            "message": f"Plugin build started for {plugin_name}",
            "task_id": task.id,
            "plugin_name": plugin_name,
            "status": "PENDING"
        }

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to start plugin build: {str(e)}"
        )


async def check_plugin_image(*, plugin_name: str) -> dict:
    """Check whether a plugin's Docker image exists locally."""
    try:
        # Get the actual plugin path (can be either official or local for checking)
        plugin_folder, source = get_plugin_path(plugin_name)

        # Official plugins don't have Docker images built locally
        if source == "official":
            return {
                "plugin_name": plugin_name,
                "image_exists": False,
                "image_tag": None,
                "source": source,
                "message": "Official plugins don't have local Docker images"
            }

        # 플러그인 폴더가 존재하는지 확인
        if not os.path.exists(plugin_folder):
            raise HTTPException(
                status_code=404,
                detail=f"Plugin folder not found: {plugin_name}"
            )

        # Dockerfile이 존재하는지 확인
        dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
        if not os.path.exists(dockerfile_path):
            raise HTTPException(
                status_code=404,
                detail=f"Dockerfile not found in plugin folder: {plugin_name}"
            )

        # Docker 이미지 존재 여부 확인
        image_exists = plugin_utils.check_plugin_docker_image(plugin_name)

        return {
            "plugin_name": plugin_name,
            "image_exists": image_exists,
            "image_tag": f"plugin-{plugin_name.lower()}" if image_exists else None,
            "source": source
        }

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to check plugin Docker image: {str(e)}"
        )


async def update_plugin_complete(*, db: Session, plugin_name: str) -> dict:
    """Synchronize a plugin's DB record from its physical metadata/dependency files."""
    try:
        # 1. 플러그인 존재 확인
        db_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin_name).first()
        if not db_plugin:
            raise HTTPException(
                status_code=404,
                detail=f"Plugin {plugin_name} not found in database"
            )

        # Check if plugin is editable
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Cannot modify plugin '{plugin_name}'. Official plugins are read-only."
            )

        # 2. 물리적 파일에서 메타데이터 읽기
        plugin_folder, source = get_plugin_path(plugin_name, "local")
        metadata_file = os.path.join(plugin_folder, "metadata.json")

        if not os.path.exists(metadata_file):
            raise HTTPException(
                status_code=404,
                detail=f"metadata.json not found for plugin {plugin_name}"
            )

        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to read metadata.json: {str(e)}"
            )

        # 3. 의존성 파일들 읽기 및 DB 업데이트 (텍스트 파일만)
        dependency_folder = os.path.join(plugin_folder, "dependency")
        dependencies_dict = {}

        if os.path.exists(dependency_folder):
            for file_name in os.listdir(dependency_folder):
                # 텍스트 의존성 파일만 DB에 저장 (.whl, .tar.gz 제외)
                if file_name.endswith(('.txt', '.yml', '.lock')):
                    file_path = os.path.join(dependency_folder, file_name)
                    try:
                        # 텍스트 파일은 내용 저장
                        with open(file_path, 'r', encoding='utf-8') as f:
                            dependencies_dict[file_name] = f.read()
                    except Exception as e:
                        print(f"Warning: Failed to read dependency file {file_name}: {str(e)}")

        # 4. DB 업데이트
        try:
            db_plugin.dependencies = dependencies_dict if dependencies_dict else None
            db_plugin.drawflow = metadata.get('drawflow', db_plugin.drawflow)
            db_plugin.rules = metadata.get('rules', db_plugin.rules)
            db_plugin.description = metadata.get('description', db_plugin.description)

            db.commit()
            db.refresh(db_plugin)

            invalidate_all_plugin_cache()
            return {
                "message": f"Plugin {plugin_name} successfully synchronized",
                "updated_fields": {
                    "dependencies": len(dependencies_dict) if dependencies_dict else 0,
                    "rules": len(metadata.get('rules', {})),
                    "drawflow_nodes": len(metadata.get('drawflow', {}).get('drawflow', {}).get('Home', {}).get('data', {})) if metadata.get('drawflow') else 0
                }
            }

        except Exception as e:
            db.rollback()
            raise HTTPException(
                status_code=500,
                detail=f"Failed to update database: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Unexpected error during plugin synchronization: {str(e)}"
        )


async def upload_text_dependencies(*, db: Session, plugin_name: str, files: List[UploadFile]) -> dict:
    """Upload text dependency files (.txt/.yml/.lock) and sync them to the DB."""
    try:
        # Check if plugin is editable
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Cannot modify plugin '{plugin_name}'. Official plugins are read-only."
            )

        # Get the actual plugin path (should be local)
        plugin_path_str, source = get_plugin_path(plugin_name, "local")
        dependency_folder = os.path.join(plugin_path_str, "dependency")

        # 1. dependency 폴더가 존재하지 않으면 생성
        if not os.path.exists(dependency_folder):
            os.makedirs(dependency_folder)

        # 2. 플러그인 DB 정보 확인
        db_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin_name).first()
        if not db_plugin:
            raise HTTPException(
                status_code=404,
                detail=f"Plugin {plugin_name} not found in database"
            )

        # 3. 업로드된 파일들 처리
        uploaded_files = []
        dependencies_dict = db_plugin.dependencies or {}

        for file in files:
            # 텍스트 의존성 파일만 허용
            if not file.filename.endswith(('.txt', '.yml', '.lock')):
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid dependency file format: {file.filename}. Only .txt, .yml, .lock files are allowed."
                )

            # 파일 내용 읽기
            try:
                content = await file.read()
                file_content = content.decode('utf-8')

                # 물리적 파일 저장
                file_path = os.path.join(dependency_folder, file.filename)
                with open(file_path, "w", encoding='utf-8') as f:
                    f.write(file_content)

                # DB용 데이터 준비
                dependencies_dict[file.filename] = file_content
                uploaded_files.append(file.filename)

                print(f"Updated dependency file: {file.filename}")

            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to process file {file.filename}: {str(e)}"
                )

        # 4. DB 업데이트
        try:
            db_plugin.dependencies = dependencies_dict
            db.commit()
            db.refresh(db_plugin)

            invalidate_all_plugin_cache()
            return {
                "message": "Text dependency files uploaded and synchronized successfully",
                "uploaded_files": uploaded_files,
                "total_dependencies": len(dependencies_dict),
                "success": True
            }

        except Exception as e:
            db.rollback()
            raise HTTPException(
                status_code=500,
                detail=f"Failed to update database: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Unexpected error during text dependency upload: {str(e)}"
        )


async def get_build_status(*, task_id: str) -> dict:
    """Return the Celery state/info for a plugin build task."""
    try:
        task_result = AsyncResult(task_id)

        result = {
            "task_id": task_id,
            "state": task_result.state,
            "current": 0,
            "total": 100,
            "info": {}
        }

        if task_result.state == 'PENDING':
            result["info"] = {"message": "Build is waiting to start..."}
        elif task_result.state == 'RUNNING':
            if task_result.info:
                result["info"] = task_result.info
            else:
                result["info"] = {"message": "Build is in progress..."}
        elif task_result.state == 'SUCCESS':
            result["current"] = 100
            result["info"] = task_result.result
        elif task_result.state == 'FAILURE':
            result["info"] = task_result.info

        return result

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get build status: {str(e)}"
        )


async def get_build_tasks(*, db: Session, user: models.User, plugin_name: Optional[str] = None) -> dict:
    """List a user's plugin_build tasks, optionally filtered by plugin name."""
    try:
        from app.task.crud import get_user_task

        # 사용자의 태스크 목록 조회 (task_type이 'plugin_build'인 것만)
        tasks = get_user_task(db, user.id)

        # plugin_build 태스크만 필터링
        build_tasks = [task for task in tasks if task.task_type == "plugin_build"]

        # plugin_name으로 필터링 (제공된 경우)
        if plugin_name:
            build_tasks = [task for task in build_tasks if task.plugin_name == plugin_name]

        # 태스크 정보를 dict로 변환
        task_list = []
        for task in build_tasks:
            # Celery 태스크 상태 조회
            task_result = AsyncResult(task.task_id)

            task_info = {
                "task_id": task.task_id,
                "plugin_name": task.plugin_name,
                "state": task_result.state,
                "start_time": task.start_time.isoformat() if task.start_time else None,
                "end_time": task.end_time.isoformat() if task.end_time else None,
                "status": task.status,
                "result": None,
                "error": None
            }

            # 태스크 결과 정보 추가
            if task_result.state == 'SUCCESS' and task_result.result:
                task_info["result"] = task_result.result
            elif task_result.state == 'FAILURE' and task_result.info:
                task_info["error"] = task_result.info
            elif task_result.state == 'RUNNING' and task_result.info:
                task_info["info"] = task_result.info

            task_list.append(task_info)

        # 시작 시간 역순으로 정렬
        task_list.sort(key=lambda x: x["start_time"] or "", reverse=True)

        return {
            "tasks": task_list,
            "count": len(task_list),
            "plugin_name": plugin_name
        }

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get build tasks: {str(e)}"
        )


async def cancel_build_task(*, task_id: str) -> dict:
    """Revoke (terminate) a running plugin build Celery task."""
    try:
        # Celery 태스크 취소
        celery_app.control.revoke(task_id, terminate=True)

        task_result = AsyncResult(task_id)

        invalidate_all_plugin_cache()
        return {
            "message": f"Build task {task_id} cancellation requested",
            "task_id": task_id,
            "state": task_result.state
        }

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to cancel build task: {str(e)}"
        )


async def get_build_logs(*, plugin_name: str) -> dict:
    """Return the last 100 lines of a local plugin's build log."""
    try:
        # Only local plugins have build logs
        if not is_plugin_editable(plugin_name):
            raise HTTPException(
                status_code=403,
                detail=f"Build logs are only available for local plugins. '{plugin_name}' is an official plugin."
            )

        # 로그 파일 경로 - use the BUILD_LOGS_DIR constant
        log_file = os.path.join(plugin_utils.BUILD_LOGS_DIR, f"{plugin_name.lower()}.log")

        if not os.path.exists(log_file):
            return {
                "plugin_name": plugin_name,
                "log_content": "Build log not found. Build may not have started yet.",
                "log_exists": False
            }

        # 로그 파일 읽기 (마지막 100줄만)
        try:
            with open(log_file, "r", encoding="utf-8") as f:
                lines = f.readlines()
                # 마지막 100줄만 가져오기
                if len(lines) > 100:
                    lines = lines[-100:]
                    log_content = f"... (showing last 100 lines) ...\n\n" + "".join(lines)
                else:
                    log_content = "".join(lines)
        except Exception as e:
            log_content = f"Error reading log file: {str(e)}"

        return {
            "plugin_name": plugin_name,
            "log_content": log_content,
            "log_exists": True,
            "log_file_path": log_file
        }

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get build logs: {str(e)}"
        )


async def get_plugin_versions(*, db: Session, plugin_name: str) -> dict:
    """Get available versions for an official plugin from the GitHub registry."""
    try:
        # Check if plugin exists and is official
        plugin = crud_plugin.get_plugin_by_name(db, plugin_name)
        if not plugin:
            raise HTTPException(
                status_code=404,
                detail=f"Plugin '{plugin_name}' not found"
            )

        if plugin.source != "official":
            raise HTTPException(
                status_code=400,
                detail=f"Version management is only available for official plugins. '{plugin_name}' is a local plugin."
            )

        # Get versions from GitHub Registry
        from app.plugin.registry import GitHubRegistryClient
        registry = GitHubRegistryClient()

        try:
            available_versions = registry.get_available_versions(plugin_name)
            current_version = plugin.version or "latest"

            return {
                "plugin_name": plugin_name,
                "current_version": current_version,
                "available_versions": available_versions,
                "total_versions": len(available_versions)
            }
        except Exception as e:
            logger.error(f"Failed to fetch versions for {plugin_name}: {e}")
            raise HTTPException(
                status_code=502,
                detail=f"Failed to fetch versions from GitHub Registry: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get plugin versions: {str(e)}"
        )


async def update_plugin_version(*, db: Session, plugin_name: str, version: str) -> dict:
    """Update an official plugin's version and (best-effort) pull its image."""
    try:
        # Check if plugin exists and is official
        plugin = crud_plugin.get_plugin_by_name(db, plugin_name)
        if not plugin:
            raise HTTPException(
                status_code=404,
                detail=f"Plugin '{plugin_name}' not found"
            )

        if plugin.source != "official":
            raise HTTPException(
                status_code=400,
                detail=f"Version updates are only available for official plugins. '{plugin_name}' is a local plugin."
            )

        # Verify version exists in registry
        from app.plugin.registry import GitHubRegistryClient
        registry = GitHubRegistryClient()

        try:
            available_versions = registry.get_available_versions(plugin_name)
            if version not in available_versions:
                raise HTTPException(
                    status_code=400,
                    detail=f"Version '{version}' not found for plugin '{plugin_name}'. Available versions: {', '.join(available_versions)}"
                )
        except Exception as e:
            logger.error(f"Failed to verify version {version} for {plugin_name}: {e}")
            raise HTTPException(
                status_code=502,
                detail=f"Failed to verify version in GitHub Registry: {str(e)}"
            )

        # Update plugin version in database
        old_version = plugin.version
        plugin.version = version

        try:
            db.commit()
            db.refresh(plugin)

            # Check if image exists locally and pull if needed
            import docker
            try:
                client = docker.from_env()
                image_uri = registry.get_image_uri(plugin_name, version)

                try:
                    client.images.get(image_uri)
                    logger.info(f"Image {image_uri} already exists locally")
                except docker.errors.ImageNotFound:
                    # Pull the new version
                    logger.info(f"Pulling new version {image_uri}...")
                    try:
                        client.images.pull(image_uri)
                        logger.info(f"Successfully pulled {image_uri}")
                    except docker.errors.APIError as e:
                        if "not found" in str(e).lower():
                            logger.warning(f"Image {image_uri} not found in registry, will fallback to local build if needed")
                        else:
                            logger.error(f"Failed to pull {image_uri}: {e}")

            except docker.errors.DockerException as e:
                logger.warning(f"Docker not available for image pull: {e}")

            invalidate_all_plugin_cache()
            return {
                "message": f"Plugin '{plugin_name}' version updated successfully",
                "plugin_name": plugin_name,
                "old_version": old_version,
                "new_version": version,
                "image_uri": registry.get_image_uri(plugin_name, version)
            }

        except Exception as e:
            db.rollback()
            logger.error(f"Failed to update plugin version in database: {e}")
            raise HTTPException(
                status_code=500,
                detail=f"Failed to update plugin version in database: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to update plugin version: {str(e)}"
        )


async def get_plugin_current_version(*, db: Session, plugin_name: str) -> dict:
    """Return current version info for a plugin (+ image URI for official plugins)."""
    try:
        plugin = crud_plugin.get_plugin_by_name(db, plugin_name)
        if not plugin:
            raise HTTPException(
                status_code=404,
                detail=f"Plugin '{plugin_name}' not found"
            )

        result = {
            "plugin_name": plugin_name,
            "current_version": plugin.version,
            "source": plugin.source,
            "is_editable": plugin.is_editable
        }

        # For official plugins, also provide image URI
        if plugin.source == "official":
            from app.plugin.registry import GitHubRegistryClient
            registry = GitHubRegistryClient()
            result["image_uri"] = registry.get_image_uri(plugin_name, plugin.version or "latest")

        return result

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to get plugin version: {str(e)}"
        )
