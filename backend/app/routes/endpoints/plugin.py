from fastapi import APIRouter, Depends, HTTPException, File, UploadFile, Form
from sqlalchemy.orm import Session
import os
from typing import List
import time
import json

from app.routes import dep
from app.database.schemas.plugin import PluginData, PluginCreate, PluginAssociate
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

        # Convert PluginData to PluginCreate
        plugin_upload_data = PluginCreate(
            name=plugin_data.plugin.name,
            description=plugin_data.plugin.description,
            author=current_user.username,  # Assuming the current user is the author
            plugin_path=f"./plugin/{plugin_data.plugin.name}/",  # Assuming this is the correct path
            dependencies=plugin_data.plugin.dependencyFiles,
            drawflow=plugin_data.drawflow,
            rules=plugin_data.rules,
        )

        #rules의 각 rule 안에 있는 script를 취합해서, 리스트로 변수 할당
        scripts = []
        for rule in plugin_data.rules:
            if rule.script:
                scripts.append(rule.script)

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

        # 중복 검사
        db_existing_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin_data.name).first()
        if db_existing_plugin:
            # 중복된 경우 업데이트
            db_plugin = crud_plugin.update_plugin(db=db, plugin=plugin_data, plugin_id=db_existing_plugin.id)
        else:
            # 중복되지 않은 경우 새로운 플러그인 생성
            db_plugin = crud_plugin.create_plugin(db=db, plugin=plugin_data)

        # 데이터베이스에 새로운 플러그인 생성
        db_plugin = crud_plugin.create_plugin(db=db, plugin=plugin_data)

        # 플러그인 폴더 생성
        if not os.path.exists(plugin_folder):
            os.makedirs(plugin_folder)
        
        # dependency 폴더 생성
        if not os.path.exists(dependency_folder):
            os.makedirs(dependency_folder)
        
        # dependency 파일 생성
        for dep in plugin_data.dependencies:
            dep_path = os.path.join(dependency_folder, dep.fileName)
            with open(dep_path, 'w') as f:
                f.write(dep.file)
        
        # Rule 객체를 사전으로 변환
        rules_dict = [rule.dict() for rule in plugin_data.rules]
        
        # metadata.json 파일 생성
        metadata = {
            "name": plugin_data.name,
            "author": plugin_data.author,
            "description": plugin_data.description,
            "drawflow": plugin_data.drawflow,
            "rules": rules_dict
        }
        metadata_path = os.path.join(plugin_folder, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

        return { "message": "Plugin data uploaded", "plugin": db_plugin }
    except Exception as e:
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
            plugin_list.append(plugin_dict)
        return { "plugins": plugin_list }
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