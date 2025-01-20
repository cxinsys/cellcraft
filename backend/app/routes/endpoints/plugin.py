from fastapi import APIRouter, Depends, HTTPException, File, UploadFile, Form
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
import os
from typing import List
import time
import json
import shutil

from app.routes import dep
from app.database.schemas.plugin import PluginData, PluginCreate, PluginUpdate, PluginAssociate
from app.database.crud import crud_plugin
from app.database import models
from app.common.utils import plugin_utils

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
            file.fileName: file.file for file in folder.files  # files 속성 사용
        }
        # 하위 폴더 처리
        for subFolder in folder.subFolders:  # subFolders 속성 사용
            folder_dict[subFolder.folderName] = parse_reference_folders([subFolder])
        result[folder.folderName] = folder_dict  # folderName 속성 사용
    return result

@router.post("/validation")
def validate_plugin(
    *,
    plugin_data: PluginData,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        # print(plugin_data)

        # Simulate validation process taking 5 seconds
        # print(plugin_data.rules)

        # Convert PluginInfo.dependencyFiles to a dictionary
        dependencies_dict = {dep.fileName: dep.file for dep in plugin_data.plugin.dependencyFiles}
        # reference_folders = {folder.folderName: {file.fileName: file.file for file in folder.files} for folder in plugin_data.plugin.referenceFolders}
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

        # rules_dict의 요소들의 isVisualization 값만 출력하기
        # print(rules_dict)

        # Convert PluginData to PluginCreate
        plugin_upload_data = PluginCreate(
            name=plugin_data.plugin.name,
            description=plugin_data.plugin.description,
            author=current_user.username,  # Assuming the current user is the author
            plugin_path=f"./plugin/{plugin_data.plugin.name}/",  # Assuming this is the correct path
            dependencies=dependencies_dict if dependencies_dict else None,
            reference_folders=reference_folders if reference_folders else None,
            drawflow=plugin_data.drawflow,
            rules=rules_dict,
        )

        # print(plugin_upload_data.rules)

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
        script_folder = os.path.join(plugin_folder, "scripts")

         # 1. 플러그인 폴더 및 의존성 폴더 생성
        plugin_utils.create_plugin_folder(plugin_folder)
        plugin_utils.create_dependency_folder(dependency_folder, plugin_data.dependencies)

        if plugin_data.reference_folders:
            print(plugin_data.reference_folders)
            plugin_utils.create_reference_folder(script_folder, plugin_data.reference_folders)
        
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
        plugin_utils.generate_snakemake_code(plugin_data.rules, plugin_folder)

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
        # 에러 발생 시 업로드된 파일 및 폴더 삭제
        if os.path.exists(plugin_folder):
            for item in os.listdir(plugin_folder):
                item_path = os.path.join(plugin_folder, item)
                if os.path.isfile(item_path):
                    os.remove(item_path)  # 파일 삭제
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)  # 디렉터리 삭제
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
        else:
            # 폴더가 이미 존재하는 경우 기존 파일 삭제
            for file in os.listdir(plugin_folder):
                file_path = os.path.join(plugin_folder, file)
                if os.path.isfile(file_path):  # 파일만 삭제
                    os.remove(file_path)
                elif os.path.isdir(file_path):  # 디렉터리 무시
                    print(f"Skipped directory: {file_path}")


        # 파일 저장
        for file in files:
            file_path = os.path.join(plugin_folder, file.filename)
            with open(file_path, "wb") as f:
                content = file.file.read()  # 파일 내용 읽기
                f.write(content)
                
        return {"message": "Scripts uploaded successfully", "scripts": [file.filename for file in files]}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

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