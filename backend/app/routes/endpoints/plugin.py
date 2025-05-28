from fastapi import APIRouter, Depends, HTTPException, File, UploadFile, Form
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
import os
from typing import List
import time
import json
import shutil
import pathlib
import uuid
import logging
from datetime import datetime
import logging

from app.routes import dep
from app.database.schemas.plugin import PluginData, PluginCreate, PluginUpdate, PluginAssociate
from app.database.crud import crud_plugin
from app.database import models
from app.common.utils import plugin_utils

# 로거 설정
logger = logging.getLogger(__name__)

router = APIRouter()

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

@router.post("/validation")
def validate_plugin(
    *,
    plugin_data: PluginData,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
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

        # Convert PluginData to PluginCreate
        plugin_upload_data = PluginCreate(
            name=plugin_data.plugin.name,
            description=plugin_data.plugin.description,
            author=current_user.username,
            plugin_path=f"./plugin/{plugin_data.plugin.name}/",
            dependencies=dependencies_dict if dependencies_dict else None,
            reference_folders=reference_folders if reference_folders else None,
            drawflow=plugin_data.drawflow,
            rules=rules_dict,
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
    
@router.post("/upload")
async def upload_plugin(
    *,
    db: Session = Depends(dep.get_db),
    plugin_data: PluginCreate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    plugin_folder = None
    backup_folder = None
    build_result = None
    try:
        # 1. 플러그인 폴더 생성
        plugin_folder = f"./plugin/{plugin_data.name}/"
        dependency_folder = os.path.join(plugin_folder, "dependency")
        script_folder = os.path.join(plugin_folder, "scripts")

        # 2. 기존 플러그인 확인 및 백업
        db_existing_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin_data.name).first()
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

        try:
            # 3. 새 폴더 생성
            plugin_utils.create_plugin_folder(plugin_folder)
            plugin_utils.create_dependency_folder(dependency_folder, plugin_data.dependencies)

            if plugin_data.reference_folders:
                plugin_utils.create_reference_folder(script_folder, plugin_data.reference_folders)
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
            plugin_utils.generate_snakemake_code(plugin_data.rules, plugin_folder, plugin_data.name)

            # 6. Dockerfile 생성
            dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
            plugin_utils.generate_plugin_dockerfile(plugin_folder, dockerfile_path)

            # 7. Docker 이미지 빌드
            build_result = plugin_utils.build_plugin_docker_image(
                plugin_path=plugin_folder,
                plugin_name=plugin_data.name,
            )

            if not build_result['success']:
                raise HTTPException(
                    status_code=500,
                    detail={
                        "message": build_result['message'],
                        "log_file": build_result['log_file']
                    }
                )

            # 8. 데이터베이스 업데이트
            if db_existing_plugin:
                db_plugin = crud_plugin.update_plugin(db=db, plugin=plugin_data, plugin_id=db_existing_plugin.id)
            else:
                db_plugin = crud_plugin.create_plugin(db=db, plugin=plugin_data)

            # 9. 트랜잭션 커밋
            db.commit()

            # 10. 성공 시 백업 삭제
            if backup_folder and os.path.exists(backup_folder):
                shutil.rmtree(backup_folder)

            return {
                "message": "플러그인 데이터 업로드 및 Docker 이미지 빌드 성공",
                "plugin": db_plugin,
                "build_info": {
                    "log_file": build_result['log_file'],
                    "image_tag": build_result['image_tag']
                }
            }

        except HTTPException as he:
            # HTTP 예외는 그대로 전달
            raise he
        except Exception as e:
            # 실패 시 롤백
            db.rollback()
            if os.path.exists(plugin_folder):
                shutil.rmtree(plugin_folder)
            if backup_folder and os.path.exists(backup_folder):
                shutil.copytree(backup_folder, plugin_folder)
                shutil.rmtree(backup_folder)
            
            error_detail = {
                "message": f"플러그인 업데이트 실패: {str(e)}",
            }
            if build_result and 'log_file' in build_result:
                error_detail["log_file"] = build_result['log_file']
            
            raise HTTPException(status_code=500, detail=error_detail)

    except HTTPException as he:
        raise he
    except Exception as e:
        if plugin_folder and os.path.exists(plugin_folder):
            shutil.rmtree(plugin_folder)
        if backup_folder and os.path.exists(backup_folder):
            shutil.rmtree(backup_folder)
        
        error_detail = {
            "message": f"플러그인 업로드 중 예기치 않은 오류: {str(e)}",
        }
        if build_result and 'log_file' in build_result:
            error_detail["log_file"] = build_result['log_file']
        
        raise HTTPException(status_code=500, detail=error_detail)
            
    
@router.post("/upload_scripts")
async def upload_scripts(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    logger = logging.getLogger(__name__)
    plugin_path = pathlib.Path(f"./plugin/{plugin_name}")

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
            
            # scripts_previous_dir 생성
            scripts_previous_dir.mkdir(parents=True, exist_ok=True)
            
            # scripts_dir의 모든 항목을 순회하면서 스크립트 파일만 이동
            for item in scripts_dir.iterdir():
                if item.is_file():  # 파일인 경우에만 이동
                    logger.info(f"Moving file {item.name} to {scripts_previous_dir}")
                    shutil.move(str(item), str(scripts_previous_dir / item.name))
                # 폴더는 그대로 유지 (reference 폴더 등)
            
            logger.info(f"Backup of script files to {scripts_previous_dir} completed.")
        else:
            logger.info(f"No current scripts directory found at {scripts_dir}. Skipping backup step.")

        # --- 핵심 트랜잭션: 스테이징 디렉토리를 최종 디렉토리로 승격 ---
        try:
            # 4. 스테이징 디렉토리를 최종 스크립트 디렉토리로 이동 (rename)
            logger.info(f"Moving staging scripts {scripts_staging_dir} to final target {scripts_dir}")
            shutil.move(str(scripts_staging_dir), str(scripts_dir))
            # 성공 시, scripts_staging_dir 경로는 더 이상 존재하지 않음
            logger.info(f"Successfully moved {scripts_staging_dir} to {scripts_dir}. Script update complete.")
            
            # `scripts_previous_dir`는 이제 바로 이전 버전의 백업이므로 삭제하지 않고 유지합니다.
            # 다음 업데이트 시, 현재의 `scripts_dir`가 새로운 `scripts_previous_dir`가 됩니다.

            return {
                "message": "Scripts uploaded successfully",
                "scripts_path": str(scripts_dir),
                "scripts": [file.filename for file in files]
            }

        except Exception as e_promote:
            logger.error(f"Failed to move staging directory {scripts_staging_dir} to {scripts_dir}: {str(e_promote)}")
            
            # 롤백 시도
            logger.info("Attempting rollback...")
            try:
                # 부분적으로 생성되었을 수 있는 새 scripts_dir 삭제
                if scripts_dir.exists(): # 실패한 새 버전일 가능성
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
                if scripts_staging_dir.exists(): # 스테이징은 확실히 정리
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
        raise he # 원래의 HTTPException을 다시 발생

    except Exception as e_unexpected:
        logger.critical(f"Unexpected error during script upload for plugin {plugin_name}: {str(e_unexpected)}", exc_info=True)
        # 예상치 못한 오류 발생 시, 최소한 스테이징 디렉토리는 정리
        if scripts_staging_dir.exists():
            try:
                shutil.rmtree(scripts_staging_dir)
                logger.info(f"Cleaned up staging directory {scripts_staging_dir} due to unexpected error.")
            except Exception as e_cleanup_staging:
                logger.error(f"Failed to cleanup staging directory {scripts_staging_dir} during unexpected error handling: {str(e_cleanup_staging)}")
        
        # scripts_dir 또는 scripts_previous_dir를 여기서 무조건 삭제하는 것은 위험할 수 있음
        # (예: 오류가 매우 초기에 발생하여 기존 운영 스크립트가 아직 안전한 상태일 경우)

        raise HTTPException(
            status_code=500,
            detail=f"An unexpected error occurred during script upload: {str(e_unexpected)}"
        )

@router.post("/upload_package")
async def upload_package(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
    ):
    try:
        # 1. 임시 폴더 생성
        temp_folder = f"./plugin/{plugin_name}/dependency_temp/"
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

        # 2. 패키지 파일을 임시 폴더에 저장
        for file in files:
            if not file.filename.endswith(('.whl', '.tar.gz')):
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid package file format: {file.filename}"
                )
            
            file_path = os.path.join(temp_folder, file.filename)
            try:
                with open(file_path, "wb") as f:
                    content = file.file.read()
                    f.write(content)
            except Exception as e:
                # 실패 시 임시 폴더 삭제
                if os.path.exists(temp_folder):
                    shutil.rmtree(temp_folder)
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to save package file {file.filename}: {str(e)}"
                )

        # 3. 기존 dependency 폴더 백업
        dependency_folder = f"./plugin/{plugin_name}/dependency/"
        backup_folder = None
        if os.path.exists(dependency_folder):
            backup_folder = f"{dependency_folder}_backup_{int(time.time())}"
            shutil.move(dependency_folder, backup_folder)

        try:
            # 4. 임시 폴더를 dependency 폴더로 이동
            shutil.move(temp_folder, dependency_folder)
            
            # 5. 성공 시 백업 폴더 삭제
            if backup_folder and os.path.exists(backup_folder):
                shutil.rmtree(backup_folder)

            return {"message": "Package uploaded successfully", "packages": [file.filename for file in files]}

        except Exception as e:
            # 실패 시 복구
            if os.path.exists(dependency_folder):
                shutil.rmtree(dependency_folder)
            if backup_folder and os.path.exists(backup_folder):
                shutil.move(backup_folder, dependency_folder)
            raise HTTPException(
                status_code=500,
                detail=f"Failed to move package files: {str(e)}"
            )

    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Unexpected error during package upload: {str(e)}"
        )

@router.get("/reference_folders/{plugin_name}")
def get_reference_folders(
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # 폴더 경로 설정
        folder_path = f"./plugin/{plugin_name}/scripts"

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

@router.get("/package/{plugin_name}")
def get_package_files(
    plugin_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # 폴더 경로 설정
        folder_path = f"./plugin/{plugin_name}/dependency"

        # 파일 리스트 중에서 .whl, .tar.gz 파일만 가져와서 반환
        package_files = [file for file in os.listdir(folder_path) if file.endswith((".whl", ".tar.gz"))]
        return {"package_files": package_files}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/file/{plugin_name}/{file_name}")
def get_plugin_file(
    plugin_name: str,
    file_name: str,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # 폴더 경로
        base_folder_path = f"./plugin/{plugin_name}"
        
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

@router.get("/list")
def list_plugins(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        plugins = crud_plugin.get_plugins(db)
        plugin_list = []
        for plugin in plugins:
            plugin_dict = plugin.__dict__
            plugin_dict['users'] = [user.__dict__ for user in plugin.users]

            # Convert rules dictionary to an array
            if 'rules' in plugin_dict:
                rules_dict = plugin_dict['rules']
                rules_array = [rules_dict[str(i)] for i in range(len(rules_dict))]
                plugin_dict['rules'] = rules_array
            
            plugin_list.append(plugin_dict)
        
        return {"plugins": plugin_list}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# plugin_name을 받아서 해당 플러그인의 정보를 반환
@router.get("/info/{plugin_name}")
def get_plugin_info(
    *,
    plugin_name: str,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # Get the plugin by name
        plugin = crud_plugin.get_plugin_by_name(db, plugin_name)
        
        if plugin is None:
            raise HTTPException(status_code=404, detail="Plugin not found")
        
        return { "plugin": plugin }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/associate")
def associate_plugin(
    *,
    db: Session = Depends(dep.get_db),
    pluginInfo: PluginAssociate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # Associate the plugin with the current user
        crud_plugin.associate_user_plugin(db, current_user.id, pluginInfo.plugin_id)

        # Get the plugin by ID
        plugin = crud_plugin.get_plugin_by_id(db, pluginInfo.plugin_id)

        return { "message": "Plugin associated with user successfully", "plugin": plugin }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@router.post("/dissociate")
def dissociate_plugin(
    *,
    db: Session = Depends(dep.get_db),
    pluginInfo: PluginAssociate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        # Dissociate the plugin with the current user
        crud_plugin.dissociate_user_plugin(db, current_user.id, pluginInfo.plugin_id)

        # Get the plugin by ID
        plugin = crud_plugin.get_plugin_by_id(db, pluginInfo.plugin_id)

        return { "message": "Plugin dissociated from user successfully", "plugin": plugin }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@router.get("/template/{plugin_id}")
def get_plugin_template(
    *,
    db: Session = Depends(dep.get_db),
    plugin_id: int,
    current_user: models.User = Depends(dep.get_current_active_user),
):
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
        
        return { "drawflow": drawflow }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))