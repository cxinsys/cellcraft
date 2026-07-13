"""Plugin file & path helpers (split from ``plugin/utils.py`` in PR-9).

Path resolution / sanitization for the official & local plugin trees, plus the
folder/file scaffolding used when a plugin is uploaded (plugin folder,
dependency folder, reference folders, metadata.json). Function signatures and
behavior are unchanged from the original ``utils.py`` implementation.
"""
import os
import re
import json
import shutil
import tempfile
from typing import List, Dict, Optional, Tuple

from fastapi import HTTPException

# Plugin directory structure constants
OFFICIAL_PLUGINS_DIR = "./plugin/official"
LOCAL_PLUGINS_DIR = "./plugin/local"


def sanitize_plugin_name(plugin_name: str) -> str:
    """
    Sanitize plugin name to prevent path traversal attacks.
    
    Args:
        plugin_name (str): Raw plugin name from user input
        
    Returns:
        str: Sanitized plugin name
        
    Raises:
        HTTPException: If plugin name contains invalid characters
    """
    if not plugin_name or not isinstance(plugin_name, str):
        raise HTTPException(status_code=400, detail="Plugin name must be a non-empty string")
    
    # Check for path traversal attempts
    if ".." in plugin_name or "/" in plugin_name or "\\" in plugin_name:
        raise HTTPException(status_code=400, detail="Plugin name contains invalid path characters")
    
    # Only allow alphanumeric, underscore, hyphen, and dots
    if not re.match(r"^[a-zA-Z0-9._-]+$", plugin_name):
        raise HTTPException(status_code=400, detail="Plugin name contains invalid characters")
    
    return plugin_name

def get_plugin_path(plugin_name: str, source: Optional[str] = None) -> Tuple[str, str]:
    """
    Get the full path to a plugin directory and determine its source.
    
    Args:
        plugin_name (str): Name of the plugin
        source (Optional[str]): If specified, look only in this source ("official" or "local")
    
    Returns:
        Tuple[str, str]: (plugin_path, actual_source)
        
    Raises:
        HTTPException: If plugin is not found or name is invalid
    """
    # Sanitize plugin name to prevent path traversal
    safe_plugin_name = sanitize_plugin_name(plugin_name)
    
    if source == "official":
        official_path = os.path.join(OFFICIAL_PLUGINS_DIR, safe_plugin_name)
        # Verify the resolved path is within the expected directory
        official_abs = os.path.abspath(official_path)
        official_base = os.path.abspath(OFFICIAL_PLUGINS_DIR)
        if not official_abs.startswith(official_base):
            raise HTTPException(status_code=400, detail="Invalid plugin path")
        if os.path.exists(official_path):
            return official_path, "official"
        raise HTTPException(status_code=404, detail=f"Official plugin '{safe_plugin_name}' not found")
    
    elif source == "local":
        local_path = os.path.join(LOCAL_PLUGINS_DIR, safe_plugin_name)
        # Verify the resolved path is within the expected directory
        local_abs = os.path.abspath(local_path)
        local_base = os.path.abspath(LOCAL_PLUGINS_DIR)
        if not local_abs.startswith(local_base):
            raise HTTPException(status_code=400, detail="Invalid plugin path")
        if os.path.exists(local_path):
            return local_path, "local"
        raise HTTPException(status_code=404, detail=f"Local plugin '{safe_plugin_name}' not found")
    
    else:
        # No source specified - prioritize local over official
        local_path = os.path.join(LOCAL_PLUGINS_DIR, safe_plugin_name)
        local_abs = os.path.abspath(local_path)
        local_base = os.path.abspath(LOCAL_PLUGINS_DIR)
        if local_abs.startswith(local_base) and os.path.exists(local_path):
            return local_path, "local"
        
        official_path = os.path.join(OFFICIAL_PLUGINS_DIR, safe_plugin_name)
        official_abs = os.path.abspath(official_path)
        official_base = os.path.abspath(OFFICIAL_PLUGINS_DIR)
        if official_abs.startswith(official_base) and os.path.exists(official_path):
            return official_path, "official"
        
        raise HTTPException(status_code=404, detail=f"Plugin '{safe_plugin_name}' not found in either local or official directories")

def list_available_plugins() -> Dict[str, List[str]]:
    """
    List all available plugins from both official and local directories.
    
    Returns:
        Dict[str, List[str]]: {"official": [...], "local": [...]}
    """
    result = {"official": [], "local": []}
    
    # List official plugins
    if os.path.exists(OFFICIAL_PLUGINS_DIR):
        try:
            official_plugins = [
                item for item in os.listdir(OFFICIAL_PLUGINS_DIR)
                if os.path.isdir(os.path.join(OFFICIAL_PLUGINS_DIR, item))
                and not item.startswith('.')
            ]
            result["official"] = sorted(official_plugins)
        except Exception as e:
            print(f"Warning: Could not list official plugins: {e}")
    
    # List local plugins
    if os.path.exists(LOCAL_PLUGINS_DIR):
        try:
            local_plugins = [
                item for item in os.listdir(LOCAL_PLUGINS_DIR)
                if os.path.isdir(os.path.join(LOCAL_PLUGINS_DIR, item))
                and not item.startswith('.')
            ]
            result["local"] = sorted(local_plugins)
        except Exception as e:
            print(f"Warning: Could not list local plugins: {e}")
    
    return result

def resolve_plugin_file_path(plugin_name: str, relative_path: str, source: Optional[str] = None) -> str:
    """
    Resolve a file path within a plugin directory.
    
    Args:
        plugin_name (str): Name of the plugin
        relative_path (str): Relative path within the plugin directory
        source (Optional[str]): Plugin source ("official" or "local")
    
    Returns:
        str: Full path to the file
        
    Raises:
        HTTPException: If plugin or file is not found
    """
    plugin_path, _ = get_plugin_path(plugin_name, source)
    full_path = os.path.join(plugin_path, relative_path)
    
    if not os.path.exists(full_path):
        raise HTTPException(
            status_code=404, 
            detail=f"File '{relative_path}' not found in plugin '{plugin_name}'"
        )
    
    return full_path

def is_plugin_editable(plugin_name: str, source: Optional[str] = None) -> bool:
    """
    Check if a plugin is editable (i.e., it's a local plugin).
    
    Args:
        plugin_name (str): Name of the plugin
        source (Optional[str]): Plugin source ("official" or "local")
    
    Returns:
        bool: True if plugin is editable, False otherwise
    """
    try:
        _, actual_source = get_plugin_path(plugin_name, source)
        return actual_source == "local"
    except HTTPException:
        return False

def ensure_local_plugins_dir():
    """
    Ensure the local plugins directory exists.
    """
    if not os.path.exists(LOCAL_PLUGINS_DIR):
        os.makedirs(LOCAL_PLUGINS_DIR, exist_ok=True)
        print(f"Created local plugins directory: {LOCAL_PLUGINS_DIR}")


def create_plugin_folder(plugin_folder: str):
    """
    Create the plugin folder if it doesn't exist.
    Preserves existing dependency folder and its contents.
    Only works for local plugins (official plugins are read-only).
    
    Parameters:
        plugin_folder (str): Path to the plugin folder.
    
    Raises:
        HTTPException: If there's an error creating the folder.
    """
    
    try:
        # Ensure this is in the local plugins directory
        ensure_local_plugins_dir()
        
        # Validate that we're not trying to modify official plugins
        if plugin_folder.startswith(OFFICIAL_PLUGINS_DIR):
            raise HTTPException(
                status_code=403,
                detail="Cannot modify official plugins. Official plugins are read-only."
            )
        
        dependency_folder = os.path.join(plugin_folder, "dependency")
        
        if os.path.exists(plugin_folder):
            # dependency 폴더가 있으면 백업
            if os.path.exists(dependency_folder):
                # 임시 디렉터리에 dependency 폴더 백업
                temp_backup_dir = tempfile.mkdtemp()
                dependency_backup = os.path.join(temp_backup_dir, "dependency_backup")
                
                try:
                    shutil.copytree(dependency_folder, dependency_backup)
                    print(f"Backed up dependency folder to: {dependency_backup}")
                    
                    # 기존 폴더 삭제
                    shutil.rmtree(plugin_folder)
                    print(f"Removed existing plugin folder: {plugin_folder}")
                    
                    # 새 플러그인 폴더 생성
                    os.makedirs(plugin_folder)
                    print(f"Created plugin folder: {plugin_folder}")
                    
                    # dependency 폴더 복원
                    shutil.copytree(dependency_backup, dependency_folder)
                    print(f"Restored dependency folder from backup")
                    
                    # 백업된 파일 목록 출력 (디버깅용)
                    restored_files = os.listdir(dependency_folder)
                    print(f"Restored dependency files: {restored_files}")
                    
                finally:
                    # 임시 백업 디렉터리 정리
                    if os.path.exists(temp_backup_dir):
                        shutil.rmtree(temp_backup_dir)
                        print(f"Cleaned up temporary backup directory")
            else:
                # dependency 폴더가 없는 경우
                shutil.rmtree(plugin_folder)
                print(f"Removed existing plugin folder: {plugin_folder}")
                os.makedirs(plugin_folder)
                print(f"Created plugin folder: {plugin_folder}")
        else:
            # 플러그인 폴더가 없는 경우 새로 생성
            os.makedirs(plugin_folder)
            print(f"Created plugin folder: {plugin_folder}")
            
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to create plugin folder: {str(e)}"
        )

def create_dependency_folder(dependency_folder: str, dependencies: dict):
    """
    Create the dependency folder and add dependency files.
    Only updates/adds files provided in dependencies dict, preserves existing package files (.whl, .tar.gz).

    Parameters:
        dependency_folder (str): Path to the dependency folder.
        dependencies (dict): A dictionary where the keys are file names and values are file contents.
    
    Raises:
        HTTPException: If there's an error creating the folder or files.
    """
    try:
        # 기존 패키지 파일들 백업
        existing_package_files = []
        if os.path.exists(dependency_folder):
            print(f"Using existing dependency folder: {dependency_folder}")
            # 기존 파일들 목록 출력 (디버깅용)
            existing_files = os.listdir(dependency_folder)
            print(f"Existing dependency files: {existing_files}")
            
            # 기존 패키지 파일들(.whl, .tar.gz) 찾기
            for file_name in existing_files:
                if file_name.endswith(('.whl', '.tar.gz')):
                    file_path = os.path.join(dependency_folder, file_name)
                    if os.path.isfile(file_path):
                        # 파일 내용을 읽어서 백업
                        with open(file_path, 'rb') as f:
                            file_content = f.read()
                        existing_package_files.append((file_name, file_content))
                        print(f"Backed up existing package file: {file_name}")
        else:
            # dependency 폴더가 없으면 생성
            os.makedirs(dependency_folder)
            print(f"Created dependency folder: {dependency_folder}")

        # dependencies 딕셔너리에 있는 파일들만 업데이트/추가
        if dependencies:
            if not isinstance(dependencies, dict):
                raise ValueError("Dependencies must be a dictionary")
            
            for file_name, file_content in dependencies.items():
                if not isinstance(file_name, str):
                    raise ValueError(f"Invalid file name: {file_name}")
                if not isinstance(file_content, str):
                    raise ValueError(f"Invalid file content for {file_name}")
                
                dep_path = os.path.join(dependency_folder, file_name)
                
                # 기존 파일이 있으면 업데이트, 없으면 새로 생성
                action = "Updated" if os.path.exists(dep_path) else "Created"
                
                with open(dep_path, 'w') as f:
                    f.write(file_content)
                print(f"{action} dependency file: {dep_path}")
        else:
            print("No dependency files to update")

        # 백업된 패키지 파일들 복원
        for file_name, file_content in existing_package_files:
            package_path = os.path.join(dependency_folder, file_name)
            if not os.path.exists(package_path):  # 새로운 의존성 파일에 포함되지 않은 경우에만 복원
                with open(package_path, 'wb') as f:
                    f.write(file_content)
                print(f"Restored existing package file: {file_name}")
            else:
                print(f"Package file already exists, skipping restore: {file_name}")

        # 최종 파일 목록 출력 (디버깅용)
        final_files = os.listdir(dependency_folder) if os.path.exists(dependency_folder) else []
        print(f"Final dependency files: {final_files}")

    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to create dependency folder or files: {str(e)}"
        )

def create_reference_folder(script_folder: str, reference_folders: dict):
    """
    Create folders and files in the reference folder structure.

    Parameters:
        script_folder (str): Base path where folders and files will be created.
        reference_folders (dict): A nested dictionary representing folder and file structures.

    Raises:
        HTTPException: If there's an error creating the folders or files.
    """
    try:
        if os.path.exists(script_folder):
            shutil.rmtree(script_folder)
        os.makedirs(script_folder)

        def process_folder(current_path: str, folder_data: dict):
            if not isinstance(folder_data, dict):
                raise ValueError(f"Invalid folder data format at {current_path}")

            for name, content in folder_data.items():
                if not isinstance(name, str):
                    raise ValueError(f"Invalid folder/file name at {current_path}")

                if name == "subFolders" and isinstance(content, list):
                    for sub_folder in content:
                        if not isinstance(sub_folder, dict):
                            raise ValueError("Invalid subfolder format")
                        for sub_name, sub_data in sub_folder.items():
                            folder_path = os.path.join(current_path, sub_name)
                            os.makedirs(folder_path, exist_ok=True)
                            process_folder(folder_path, sub_data)
                elif isinstance(content, dict):
                    folder_path = os.path.join(current_path, name)
                    os.makedirs(folder_path, exist_ok=True)
                    process_folder(folder_path, content)
                else:
                    if not isinstance(content, str):
                        raise ValueError(f"Invalid file content for {name}")
                    file_path = os.path.join(current_path, name)
                    with open(file_path, 'w') as f:
                        f.write(content)

        process_folder(script_folder, reference_folders)

    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to create reference folder structure: {str(e)}"
        )

def get_reference_folders_list(folder_path: str) -> list:
    """
    Get a list of folder names in the specified folder path.

    Parameters:
        folder_path (str): Path to the folder.

    Returns:
        list: A list of folder names in the specified folder path.
    """
    if not os.path.exists(folder_path):
        raise HTTPException(status_code=404, detail="Folder not found")

    return [item for item in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, item))]


def get_reference_folder(folder_path: str) -> dict:
    """
    Recursively retrieve the structure of a single folder.

    Parameters:
        folder_path (str): Path to the folder.

    Returns:
        dict: Folder structure.
    """
    if not os.path.exists(folder_path):
        raise HTTPException(status_code=404, detail="Folder not found")

    folder_structure = {
        "folderName": os.path.basename(folder_path),
        "files": [],
        "subFolders": []
    }

    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path):
            folder_structure["files"].append({
                "name": item,
                "type": ""
            })
        elif os.path.isdir(item_path):
            folder_structure["subFolders"].append(get_reference_folder(item_path))

    return folder_structure

def create_metadata_file(plugin_folder: str, metadata: dict):
    """
    Create the metadata.json file in the plugin folder.
    
    Parameters:
        plugin_folder (str): Path to the plugin folder.
        metadata (dict): Metadata dictionary to be saved as metadata.json.
    """
    metadata_path = os.path.join(plugin_folder, "metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=4)
    print(f"Metadata file created at {metadata_path}")
