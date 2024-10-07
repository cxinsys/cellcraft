from fastapi import APIRouter, Depends, HTTPException, File, UploadFile, Form
from sqlalchemy.orm import Session
import os
from typing import List
import time
import json

from app.routes import dep
from app.database.schemas.plugin import PluginData, PluginCreate, PluginUpdate, PluginAssociate
from app.database.crud import crud_plugin
from app.database import models
from app.common.utils import plugin_utils

router = APIRouter()

@router.post("/validation")
def validate_plugin(
    *,
    plugin_data: PluginData,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        print(plugin_data)

        # Simulate validation process taking 5 seconds
        time.sleep(5)

        # Convert PluginInfo.dependencyFiles to a dictionary
        dependencies_dict = {dep.fileName: dep.file for dep in plugin_data.plugin.dependencyFiles}
        rules_dict = {index: rule for index, rule in enumerate(plugin_data.rules)}

        # Convert PluginData to PluginCreate
        plugin_upload_data = PluginCreate(
            name=plugin_data.plugin.name,
            description=plugin_data.plugin.description,
            author=current_user.username,  # Assuming the current user is the author
            plugin_path=f"./plugin/{plugin_data.plugin.name}/",  # Assuming this is the correct path
            dependencies=dependencies_dict if dependencies_dict else None,
            drawflow=plugin_data.drawflow,
            rules=rules_dict,
        )

        # Collect scripts from rules
        scripts = [rule.script for rule in plugin_data.rules if rule.script]

        return { "message": "Plugin data validated", "plugin": plugin_upload_data, "scripts": scripts }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@router.post("/upload")
def upload_plugin(
    *,
    db: Session = Depends(dep.get_db),
    plugin_data: PluginCreate,
    current_user: models.User = Depends(dep.get_current_active_user),
):
    try:
        plugin_folder = f"./plugin/{plugin_data.name}/"
        dependency_folder = os.path.join(plugin_folder, "dependency")

         # 1. 플러그인 폴더 및 의존성 폴더 생성
        plugin_utils.create_plugin_folder(plugin_folder)
        plugin_utils.create_dependency_folder(dependency_folder, plugin_data.dependencies)
        
        # 2. 메타데이터 생성
        metadata = {
            "name": plugin_data.name,
            "author": plugin_data.author,
            "description": plugin_data.description,
            "drawflow": plugin_data.drawflow,
            "rules": plugin_data.rules
        }
        plugin_utils.create_metadata_file(plugin_folder, metadata)
        
        # 3. Snakefile 생성
        snakefile_path = os.path.join(plugin_folder, "Snakefile")
        plugin_utils.generate_snakemake_code(plugin_data.rules, snakefile_path)

        # 4. 의존성 설치
        # requirements.txt, environment.yml, renv.lock 파일을 확인 후 설치
        # for dependency_file in ['requirements.txt', 'environment.yml', 'environment.yaml', 'renv.lock']:
        #     dependency_file_path = os.path.join(dependency_folder, dependency_file)
        #     if os.path.exists(dependency_file_path):
        #         print(f"Installing dependencies from {dependency_file}...")
        #         plugin_utils.install_dependencies(dependency_file_path)

        # 5. 중복 검사 및 플러그인 생성 또는 업데이트
        db_existing_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin_data.name).first()
        if db_existing_plugin:
            # 중복된 경우 업데이트
            db_plugin = crud_plugin.update_plugin(db=db, plugin=plugin_data, plugin_id=db_existing_plugin.id)
        else:
            # 중복되지 않은 경우 새로운 플러그인 생성
            db_plugin = crud_plugin.create_plugin(db=db, plugin=plugin_data)

        return { "message": "Plugin data uploaded", "plugin": db_plugin }
    except Exception as e:
        # 에러 발생 시 업로드된 파일이 포함된 플러그인 폴더 삭제
        if os.path.exists(plugin_folder):
            os.rmdir(plugin_folder)

        raise HTTPException(status_code=500, detail=str(e))
    
@router.post("/upload_scripts")
def upload_scripts(
    plugin_name: str = Form(...),
    files: List[UploadFile] = File(...),
    current_user: models.User = Depends(dep.get_current_active_user),
    ):
    try:
        print(plugin_name)
        # 스크립트 파일을 저장할 폴더 경로
        plugin_folder = f"./plugin/{plugin_name}/scripts/"
        
        # 폴더가 없으면 생성
        if not os.path.exists(plugin_folder):
            os.makedirs(plugin_folder)

        # 파일 저장
        for file in files:
            file_path = os.path.join(plugin_folder, file.filename)
            with open(file_path, "wb") as f:
                content = file.file.read()  # 파일 내용 읽기
                f.write(content)
                
        return {"message": "Scripts uploaded successfully", "scripts": [file.filename for file in files]}
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