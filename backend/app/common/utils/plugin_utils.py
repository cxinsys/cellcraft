import json
import os
import re
import subprocess
import yaml
import shutil
import numpy as np
import time
import platform
from fastapi import HTTPException
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import docker
from datetime import datetime

# Plugin directory structure constants
OFFICIAL_PLUGINS_DIR = "./plugin/official"
LOCAL_PLUGINS_DIR = "./plugin/local"
BUILD_LOGS_DIR = "./plugin/build_logs"

# Platform detection and Docker build utilities
def detect_target_platform():
    """
    시스템 아키텍처를 자동 감지하여 Docker 플랫폼 문자열 반환
    
    Returns:
        str: Docker 플랫폼 문자열 (예: "linux/amd64", "linux/arm64")
    """
    machine = platform.machine().lower()
    system = platform.system().lower()
    
    # 아키텍처 매핑 테이블
    arch_mapping = {
        'x86_64': 'amd64',
        'amd64': 'amd64',
        'aarch64': 'arm64',
        'arm64': 'arm64',
        'armv7l': 'arm/v7',
        'armv6l': 'arm/v6'
    }
    
    # 시스템이 Linux가 아닌 경우 기본값 반환
    if system != 'linux':
        system = 'linux'  # Docker 컨테이너는 Linux 기반
    
    arch = arch_mapping.get(machine, 'amd64')  # 알 수 없는 아키텍처는 amd64 기본값
    return f"{system}/{arch}"

def detect_platform_via_docker(client):
    """
    Docker 데몬에서 플랫폼 정보 추출
    
    Args:
        client: Docker 클라이언트 객체
        
    Returns:
        str: Docker 플랫폼 문자열, 실패 시 "linux/amd64"
    """
    try:
        info = client.info()
        architecture = info.get('Architecture', 'x86_64')
        os_type = info.get('OSType', 'linux').lower()
        
        arch_map = {
            'x86_64': 'amd64',
            'aarch64': 'arm64',
            'arm64': 'arm64'
        }
        
        docker_arch = arch_map.get(architecture, 'amd64')
        return f"{os_type}/{docker_arch}"
    except Exception:
        return "linux/amd64"  # fallback

def get_target_platform(client=None):
    """
    우선순위에 따른 플랫폼 감지
    
    우선순위:
    1. 환경변수 DOCKER_DEFAULT_PLATFORM
    2. Docker 데몬 정보 
    3. Python 시스템 정보
    
    Args:
        client: Docker 클라이언트 객체 (선택사항)
        
    Returns:
        str: 감지된 Docker 플랫폼 문자열
    """
    # 1순위: 환경변수 직접 지정
    if env_platform := os.environ.get('DOCKER_DEFAULT_PLATFORM'):
        return env_platform
    
    # 2순위: Docker 데몬 정보
    if client:
        docker_platform = detect_platform_via_docker(client)
        if docker_platform != "linux/amd64":  # 기본값이 아닌 경우 우선 사용
            return docker_platform
    
    # 3순위: Python 시스템 정보
    return detect_target_platform()

def check_gpu_requirement(plugin_path: str) -> bool:
    """
    플러그인의 GPU 요구사항 확인
    
    Args:
        plugin_path (str): 플러그인 폴더 경로
        
    Returns:
        bool: GPU 필요 여부
    """
    try:
        metadata_path = os.path.join(plugin_path, "metadata.json")
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                return metadata.get('requiresGPU', False)
    except Exception:
        pass
    
    # metadata에서 확인 불가능한 경우, 기본값 False (CPU 전용)
    return False

def get_optimal_platform(plugin_path: str, client=None):
    """
    GPU 사용 여부에 따른 최적 플랫폼 선택
    
    Args:
        plugin_path (str): 플러그인 폴더 경로
        client: Docker 클라이언트 객체 (선택사항)
        
    Returns:
        str: 최적화된 Docker 플랫폼 문자열
    """
    use_gpu = check_gpu_requirement(plugin_path)
    
    if use_gpu:
        # GPU는 현재 linux/amd64만 지원 (NVIDIA CUDA)
        return "linux/amd64"
    
    # CPU 전용은 자동 감지
    return get_target_platform(client)

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

def ensure_build_logs_dir():
    """
    Ensure the build logs directory exists.
    """
    if not os.path.exists(BUILD_LOGS_DIR):
        os.makedirs(BUILD_LOGS_DIR, exist_ok=True)
        print(f"Created build logs directory: {BUILD_LOGS_DIR}")

def generate_merged_plugin_drawflow( drawflow_data_list: List[Dict[str, Any]], plugin_name: str):
    merged_drawflow = {"drawflow": {"Home": {"data": {}}}}
    node_id = 1  # ID 초기화
    inputfile_nodes = {}  # 이미 생성된 InputFile 노드를 저장 (key: 파일명, value: 노드 ID)

    pos_x_start = 50
    pos_y_start = 50
    pos_y_increment = 100
    pos_x_increment = 200

    pos_x_datatable = pos_x_start + pos_x_increment
    pos_x_algorithm = pos_x_datatable + pos_x_increment

    pos_y = pos_y_start

    new_nodes = []
    algorithm_inputs = []

    for drawflow_data in drawflow_data_list:
        original_data = drawflow_data["drawflow"]["Home"]["data"]

        for key, node in original_data.items():
            if node["name"] == "InputFile":
                input_file = node["data"]["title"]

                # 중복된 InputFile 노드가 있으면 기존 노드 사용
                if input_file in inputfile_nodes:
                    inputfile_node_id = inputfile_nodes[input_file]
                else:
                    # 새로운 InputFile 노드 생성
                    inputfile_node_id = node_id
                    inputfile_nodes[input_file] = node_id

                    inputfile_node = {
                        "id": inputfile_node_id,
                        "name": "InputFile",
                        "data": {"title": input_file},
                        "class": "InputFile",
                        "html": "InputFile",
                        "typenode": "vue",
                        "inputs": {},
                        "outputs": {"output_1": {"connections": []}},
                        "pos_x": pos_x_start,
                        "pos_y": pos_y
                    }
                    merged_drawflow["drawflow"]["Home"]["data"][str(node_id)] = inputfile_node
                    node_id += 1
                    pos_y += pos_y_increment

                # DataTable 노드 추가
                datatable_node_id = node_id
                datatable_node = {
                    "id": datatable_node_id,
                    "name": "DataTable",
                    "data": {},
                    "class": "DataTable",
                    "html": "DataTable",
                    "typenode": "vue",
                    "inputs": {"input_1": {"connections": [{"node": str(inputfile_node_id), "input": "output_1"}]}},
                    "outputs": {"output_1": {"connections": []}},
                    "pos_x": pos_x_datatable,
                    "pos_y": pos_y_start
                }
                merged_drawflow["drawflow"]["Home"]["data"][str(node_id)] = datatable_node
                node_id += 1

                algorithm_inputs.append({"node": str(datatable_node_id), "input": "output_1"})

    # Algorithm 노드 생성
    algorithm_node = {
        "id": node_id,
        "name": "Algorithm",
        "data": {"title": plugin_name},
        "class": "Algorithm",
        "html": "Algorithm",
        "typenode": "vue",
        "inputs": {"input_1": {"connections": algorithm_inputs}},
        "outputs": {"output_1": {"connections": []}},
        "pos_x": pos_x_algorithm,
        "pos_y": pos_y_start
    }
    merged_drawflow["drawflow"]["Home"]["data"][str(node_id)] = algorithm_node

    # 생성된 drawflow 출력
    print(json.dumps(merged_drawflow, indent=2))
    return merged_drawflow

def generate_plugin_drawflow_template(drawflow_data: Dict[str, Any], plugin_name: str):
    original_data = drawflow_data['drawflow']['Home']['data']

    # 각 노드의 ID를 생성합니다.
    node_id = max(int(key) for key in original_data.keys()) + 1

    # 새로운 drawflow 구조를 초기화합니다.
    new_drawflow = {
        "drawflow": {
            "Home": {
                "data": {}
            }
        }
    }

    # 노드별 위치 초기값 설정
    pos_x_start = 50
    pos_y_start = 50
    pos_y_increment = 100
    pos_x_increment = 200

    pos_x_inputfile = pos_x_start
    pos_x_datatable = pos_x_start + pos_x_increment
    pos_x_scatterplot = pos_x_datatable + pos_x_increment
    pos_x_algorithm = pos_x_scatterplot + pos_x_increment
    pos_x_resultfile = pos_x_algorithm + pos_x_increment
    pos_x_visualization = pos_x_resultfile + pos_x_increment

    pos_y_inputfile = pos_y_start
    pos_y_datatable = pos_y_start
    pos_y_scatterplot = pos_y_start
    pos_y_algorithm = pos_y_start
    pos_y_resultfile = pos_y_start
    pos_y_visualization = pos_y_start

    # 새로운 노드를 저장할 리스트 초기화
    new_nodes = []
    resultfile_nodes = []
    visualization_nodes = []

    # 기존 노드를 순회하며 새로운 노드를 생성합니다.
    for key, node in original_data.items():
        # inputs를 확인하여 connections가 비어있는 경우 새로운 노드 생성
        for input_key, input_val in node['inputs'].items():
            if node['data'].get('isVisualization', False):
                continue

            if not input_val['connections']:
                index = int(input_key.split('_')[1]) - 1
                input_file = node['data']['inputs'][index]

                # input_file과 같은 값으로 이미 생성된 InputFile 노드가 있는지 확인 후, 존재하면 continue
                if any(str(node_id_) in new_drawflow["drawflow"]["Home"]["data"] and 
                       new_drawflow["drawflow"]["Home"]["data"][str(node_id_)]["name"] == "InputFile" and 
                       new_drawflow["drawflow"]["Home"]["data"][str(node_id_)]["data"]["title"] == input_file 
                       for node_id_ in new_drawflow["drawflow"]["Home"]["data"]):
                    continue

                if input_file.endswith('.h5ad'):
                    # InputFile 노드 생성
                    inputfile_node = {
                        "id": node_id,
                        "name": "InputFile",
                        "data": {"title": input_file},
                        "class": "InputFile",
                        "html": "InputFile",
                        "typenode": "vue",
                        "inputs": {},
                        "outputs": {
                            "output_1": {
                                "connections": [{"node": str(node_id + 1), "input": "input_1"}]
                            }
                        },
                        "pos_x": pos_x_inputfile,
                        "pos_y": pos_y_inputfile
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = inputfile_node
                    inputfile_node_id = node_id
                    node_id += 1
                    pos_y_inputfile += pos_y_increment

                    # DataTable 노드 생성
                    datatable_node = {
                        "id": node_id,
                        "name": "DataTable",
                        "data": {},
                        "class": "DataTable",
                        "html": "DataTable",
                        "typenode": "vue",
                        "inputs": {
                            "input_1": {
                                "connections": [{"node": str(inputfile_node_id), "input": "output_1"}]
                            }
                        },
                        "outputs": {
                            "output_1": {
                                "connections": [{"node": str(node_id + 1), "input": "input_1"}]
                            }
                        },
                        "pos_x": pos_x_datatable,
                        "pos_y": pos_y_datatable
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = datatable_node
                    datatable_node_id = node_id
                    node_id += 1
                    pos_y_datatable += pos_y_increment

                    # ScatterPlot 노드 생성
                    scatterplot_node = {
                        "id": node_id,
                        "name": "ScatterPlot",
                        "data": {},
                        "class": "ScatterPlot",
                        "html": "ScatterPlot",
                        "typenode": "vue",
                        "inputs": {
                            "input_1": {
                                "connections": [{"node": str(datatable_node_id), "input": "output_1"}]
                            }
                        },
                        "outputs": {
                            "output_1": {
                                "connections": []
                            }
                        },
                        "pos_x": pos_x_scatterplot,
                        "pos_y": pos_y_scatterplot
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = scatterplot_node
                    new_nodes.append(scatterplot_node)
                    node_id += 1
                    pos_y_scatterplot += pos_y_increment
                else:
                    inputfile_node = {
                        "id": node_id,
                        "name": "InputFile",
                        "data": {"title": input_file},
                        "class": "InputFile",
                        "html": "InputFile",
                        "typenode": "vue",
                        "inputs": {},
                        "outputs": {
                            "output_1": {
                                "connections": []
                            }
                        },
                        "pos_x": pos_x_inputfile,
                        "pos_y": pos_y_inputfile
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = inputfile_node
                    new_nodes.append(inputfile_node)
                    node_id += 1
                    pos_y_inputfile += pos_y_increment

        # outputs를 확인하여 connections가 비어있는 경우 새로운 노드 생성
        for output_key, output_val in node['outputs'].items():
            # connections가 비어있거나, 모든 연결이 isVisualization이 True인 노드에만 연결되어 있는지 확인
            should_create_resultfile = False
            if not output_val['connections']:
                should_create_resultfile = True
            else:
                # 모든 연결이 isVisualization이 True인 노드인지 확인
                all_visualization = True
                for conn in output_val['connections']:
                    connected_node = original_data.get(conn['node'])
                    if connected_node and not connected_node.get('data', {}).get('isVisualization', False):
                        all_visualization = False
                        break
                should_create_resultfile = all_visualization

            if should_create_resultfile:
                index = int(output_key.split('_')[1]) - 1
                output_file = node['data']['outputs'][index]
                if node['data'].get('isVisualization', False):
                    visualization_node = {
                        "id": node_id,
                        "name": "Visualization",
                        "data": {"title": output_file},
                        "class": "Visualization",
                        "html": "Visualization",
                        "typenode": "vue",
                        "inputs": {
                            "input_1": {
                                "connections": []
                            }
                        },
                        "outputs": {},
                        "pos_x": pos_x_visualization,
                        "pos_y": pos_y_visualization
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = visualization_node
                    visualization_nodes.append(visualization_node)
                    node_id += 1
                    pos_y_visualization += pos_y_increment
                else:
                    resultfile_node = {
                        "id": node_id,
                        "name": "ResultFile",
                        "data": {"title": output_file},
                        "class": "ResultFile",
                        "html": "ResultFile",
                        "typenode": "vue",
                        "inputs": {
                            "input_1": {
                                "connections": []
                            }
                        },
                        "outputs": {
                            "output_1": {
                                "connections": []
                            }
                        },
                        "pos_x": pos_x_resultfile,
                        "pos_y": pos_y_resultfile
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = resultfile_node
                    resultfile_nodes.append(resultfile_node)
                    node_id += 1
                    pos_y_resultfile += pos_y_increment

    # Algorithm 노드 생성
    algorithm_node = {
        "id": node_id,
        "name": "Algorithm",
        "data": {"title": plugin_name},
        "class": "Algorithm",
        "html": "Algorithm",
        "typenode": "vue",
        "inputs": {
            "input_1": {
                "connections": []
            }
        },
        "outputs": {
            "output_1": {
                "connections": []
            }
        },
        "pos_x": pos_x_algorithm,
        "pos_y": pos_y_algorithm
    }

    for new_node in new_nodes:
        if new_node['name'] in ["InputFile", "DataTable", "ScatterPlot"]:
            algorithm_node["inputs"]["input_1"]["connections"].append({
                "node": str(new_node['id']),
                "input": "output_1"
            })
            new_node["outputs"]["output_1"]["connections"].append({
                "node": str(algorithm_node['id']),
                "input": "input_1"
            })

    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = algorithm_node
    algorithm_node_id = node_id
    node_id += 1

    # Algorithm 노드의 output과 모든 ResultFile 노드의 input을 연결
    for resultfile_node in resultfile_nodes:
        resultfile_node["inputs"]["input_1"]["connections"].append({
            "node": str(algorithm_node_id),
            "input": "output_1"
        })
        algorithm_node["outputs"]["output_1"]["connections"].append({
            "node": str(resultfile_node['id']),
            "output": "input_1"
        })

    # visualization 노드의 입력을 모든 ResultFile 노드의 출력에 연결
    for visualization_node in visualization_nodes:
        visualization_node["inputs"]["input_1"]["connections"] = [
            {"node": str(resultfile_node['id']), "input": "output_1"}
            for resultfile_node in resultfile_nodes
        ]
        for resultfile_node in resultfile_nodes:
            resultfile_node["outputs"]["output_1"]["connections"].append({
                "node": str(visualization_node['id']),
                "input": "input_1"
            })

    # 생성된 new_drawflow 데이터를 출력합니다.
    print(json.dumps(new_drawflow, indent=2))
    return new_drawflow

def normalize_param_name(name: str) -> str:
    """파라미터 이름에서 특수문자를 제거하고 언더스코어로 변환합니다."""
    # 공백을 언더스코어로 변환
    name = name.replace(' ', '_')
    # 알파벳, 숫자, 언더스코어만 허용하고 나머지는 제거
    name = re.sub(r'[^a-zA-Z0-9_]', '', name)
    return name

def normalize_string(text: str) -> str:
    """문자열에서 특수문자를 제거합니다."""
    # 알파벳, 숫자, 공백, 점(.), 언더스코어(_), 하이픈(-) 만 허용
    return re.sub(r'[^a-zA-Z0-9\s._-]', '', text)

def generate_snakemake_code(rules_data, output_folder_path, plugin_type=None, workflow_id=None, algo_id=None, viz_id=None):
    """
    Generate unified Snakemake code for both analysis and visualization plugins.
    
    Parameters:
        rules_data (dict): Rules data
        output_folder_path (str): Output folder path
        plugin_type (str): Plugin type ('ANALYSIS' or 'VISUALIZATION')
        workflow_id (str): Workflow ID for visualization plugins
        algo_id (str): Algorithm ID for input file resolution
        viz_id (str): Visualization node ID for output path
    
    Raises:
        HTTPException: If there's an error generating the Snakemake code
    """
    try:
        if not rules_data:
            raise ValueError("No rules data provided")

        snakemake_path = os.path.join(output_folder_path, "Snakefile")
        snakemake_code = "import os\n\n"  # Add import statement at the top
        
        # Configure paths based on plugin type
        # When workflow_id, algo_id, viz_id are provided, use them; otherwise use placeholders
        if plugin_type == 'VISUALIZATION':
            # For visualization plugins: input from algorithm results, output to visualization directory
            if workflow_id and algo_id and viz_id:
                # Specific values provided - for runtime execution
                input_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                output_path = f"user/{{user_name}}/workflow_{workflow_id}/visualization_{viz_id}/results"
                logs_path = f"user/{{user_name}}/workflow_{workflow_id}/visualization_{viz_id}/logs"
            else:
                # No specific values - generate template with placeholders
                input_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                output_path = "user/{user_name}/workflow_{workflow_id}/visualization_{visualization_id}/results"
                logs_path = "user/{user_name}/workflow_{workflow_id}/visualization_{visualization_id}/logs"
            unique_input_path = input_path  # Visualization plugins don't use unique inputs from data folder
        else:
            # For analysis plugins: input from data, output to algorithm results
            if workflow_id and algo_id:
                # Specific values provided - for runtime execution
                input_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                output_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                logs_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/logs"
            else:
                # No specific values - generate template with placeholders
                input_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                output_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                logs_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/logs"
            unique_input_path = "user/{user_name}/data"

        # Find all outputs across all rules to determine unique inputs not present in outputs
        all_outputs = {out for rule in rules_data.values() for out in rule['output']}
        all_inputs = {inp for rule in rules_data.values() for inp in rule['input']}
        unique_inputs = [inp for inp in all_inputs if inp not in all_outputs]
        scripts_log_word = "> {log.stdout} 2> {log.stderr}"

        # Iterate through each rule in the dictionary
        for rule_id, rule in rules_data.items():
            # All rules go into the same unified Snakefile
            print(f"Adding rule to Snakefile: {rule['name']} (plugin_type: {plugin_type})")
            
            # Start rule block with rule name
            snakemake_code += f"rule {normalize_param_name(rule['name'])}:\n"

            # Input section with optional input handling
            if 'input' in rule and rule['input']:
                def get_input_param_name(inp, rule_params):
                    for param in rule_params:
                        if param['type'] == 'inputFile' and param.get('defaultValue') == inp:
                            return normalize_param_name(param['name'])
                    return normalize_param_name(os.path.splitext(os.path.basename(inp))[0])

                input_files = ",\n        ".join([
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_path}/{{{inp}}}\"" 
                    if 'target' in get_input_param_name(inp, rule.get('parameters', [])) else
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{unique_input_path}/{{{inp}}}\"" 
                    if inp in unique_inputs else
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_path}/{inp}\""
                    for inp in rule['input']
                    if not any(param['type'] == 'optionalInputFile' and param.get('defaultValue') == inp 
                              for param in rule.get('parameters', []))
                ])
                if input_files:  # Only add input section if there are non-optional input files
                    snakemake_code += f"    input:\n        {input_files}\n"

            # Output section
            if 'output' in rule and rule['output']:
                def get_output_param_name(out, rule_params):
                    for param in rule_params:
                        if param['type'] == 'outputFile' and param.get('defaultValue') == out:
                            return normalize_param_name(param['name'])
                    return normalize_param_name(os.path.splitext(os.path.basename(out))[0])

                # Use appropriate output path based on plugin type
                output_files = ",\n        ".join([
                    f"{get_output_param_name(out, rule.get('parameters', []))}=\"{output_path}/{out}\""
                    for out in rule['output']
                ])
                snakemake_code += f"    output:\n        {output_files}\n"

            # Params section
            if 'parameters' in rule and rule['parameters']:
                param_list = []
                for param in rule['parameters']:
                    normalized_name = normalize_param_name(param['name'])
                    if param['name'] == "clusters" and param['type'] == "h5adParameter":
                        param_list.append(f'clusters=lambda wc: ";".join({{{param["name"]}}})')
                    elif param['name'] == "ScatterPlot" and param['type'] == 'string':
                        param_list.append(
                            'ScatterPlot=lambda wildcards: {ScatterPlot} if os.path.exists({ScatterPlot}) else "None"'
                        )
                    elif "UMAP" in param['name'] and param['type'] == 'h5adParameter':
                        param_list.append(
                            'UMAP_lasso=lambda wildcards: {UMAP lasso} if os.path.exists({UMAP lasso}) else "None"'
                        )
                    elif param['type'] == 'optionalInputFile':
                        param_list.append(
                            f"{normalized_name}=lambda wildcards: \"{unique_input_path}/{{{param['defaultValue']}}}\" "
                            f"if os.path.exists(\"{unique_input_path}/{{{param['defaultValue']}}}\") else \"None\""
                        )
                    elif param['type'] != 'inputFile' and param['type'] != 'outputFile':
                        param_list.append(f"{normalized_name}={{{param['name']}}}")

                if param_list:
                    param_list_str = ",\n        ".join(param_list)
                    snakemake_code += f"    params:\n        {param_list_str}\n"

            # Log section
            snakemake_code += f"    log:\n"
            snakemake_code += f"        stdout=\"{logs_path}/{normalize_param_name(rule['name'])}.stdout\",\n"
            snakemake_code += f"        stderr=\"{logs_path}/{normalize_param_name(rule['name'])}.stderr\"\n"

            # Shell section based on script type
            if 'script' in rule and rule['script']:
                # Script path 설정 (Docker 내부 경로)
                script_path = f"/scripts/{normalize_string(rule['script'])}"

                # Shell command 설정 (Python/R 분기)
                if rule['script'].endswith('.py'):
                    shell_command = f"/opt/micromamba/envs/plugin_env/bin/python {script_path}"
                elif rule['script'].endswith('.R'):
                    # Simplified R script execution using micromamba environment
                    shell_command = f"export R_HOME=/opt/micromamba/envs/plugin_env/lib/R && " \
                                  f"export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE && " \
                                  f"/opt/micromamba/envs/plugin_env/bin/Rscript {script_path}"
                else:
                    shell_command = script_path

                # Add parameters to shell command with input/output/params distinction
                if 'parameters' in rule and rule['parameters']:
                    param_list = []
                    for param in rule['parameters']:
                        normalized_name = normalize_param_name(param['name'])
                        if param['type'] == 'inputFile':
                            param_list.append(f"{{input.{normalized_name}}}")
                        elif param['type'] == 'optionalInputFile':
                            param_list.append(f"{{params.{normalized_name}}}")
                        elif param['type'] == 'outputFile':
                            param_list.append(f"{{output.{normalized_name}}}")
                        elif param['name'] == "clusters" and param['type'] == "h5adParameter":
                            param_list.append(f"'{{{normalize_string(f'params.{normalized_name}')}}}'")
                        else:
                            param_list.append(f"{{{normalize_string(f'params.{normalized_name}')}}}")

                    param_list_str = " ".join(param_list)
                    shell_command = f"{shell_command} {param_list_str}"

                snakemake_code += f"    shell:\n        \"{shell_command} {scripts_log_word}\"\n"

            snakemake_code += "\n"  # Add a newline for separation between rules

        # Write the unified Snakefile
        with open(snakemake_path, 'w') as file:
            file.write(snakemake_code)
            
        print(f"Generated unified Snakefile at: {snakemake_path} (plugin_type: {plugin_type})")

    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to generate Snakemake code: {str(e)}"
        )

def normalize_pkg_name(name: str) -> str:
    """Normalize package names for consistent comparison."""
    return name.lower().replace("-", "_")

def check_requirements_txt(requirements_file: str):
    """Check if dependencies in requirements.txt are installed."""
    
    # Check if the requirements file exists
    if not os.path.exists(requirements_file):
        raise FileNotFoundError(f"{requirements_file} not found.")
    
    # Get installed packages via pip freeze
    try:
        installed_packages = subprocess.run(
            ["/opt/conda/envs/snakemake/bin/pip", "freeze"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        ).stdout.splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Error getting installed packages: {e}")
        return {
            "installed_status": False,
            "message": "Failed to get installed packages"
        }

    # Read the requirements.txt file
    with open(requirements_file, 'r') as f:
        required_packages = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    # Parse installed packages to a dict
    installed_dict = {}
    for pkg in installed_packages:
        if '==' in pkg:
            name, version = pkg.split('==', 1)
            installed_dict[normalize_pkg_name(name)] = version

    missing = []
    # Check each requirement
    for req in required_packages:
        # Skip empty lines and comments
        if not req or req.startswith('#'):
            continue
            
        # Parse requirement line
        pkg_name = req.split('==')[0].split('>=')[0].split('<=')[0].split('~=')[0].strip()
        pkg_name = normalize_pkg_name(pkg_name)
        
        if pkg_name not in installed_dict:
            missing.append(pkg_name)

    # Return result in dict form
    if missing:
        return {
            "installed_status": False,
            "missing_packages": missing
        }
    
    return {
        "installed_status": True
    }

def check_environment_yml(environment_file: str):
    """Check if dependencies in environment.yml are installed, including pip packages."""
    
    # Check if the environment file exists
    if not os.path.exists(environment_file):
        raise FileNotFoundError(f"{environment_file} not found.")
    
    try:
        # Get installed conda packages
        conda_result = subprocess.run(
            ["/opt/conda/envs/snakemake/bin/conda", "list", "--json"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        installed_conda_packages = json.loads(conda_result.stdout)
        
        # Get installed pip packages
        pip_result = subprocess.run(
            ["/opt/conda/envs/snakemake/bin/pip", "list", "--format=json"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        installed_pip_packages = json.loads(pip_result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error getting installed packages: {e}")
        return {
            "installed_status": False,
            "message": "Failed to get installed packages"
        }
    except json.JSONDecodeError as e:
        print(f"Error parsing package list: {e}")
        return {
            "installed_status": False,
            "message": "Failed to parse package list"
        }

    # Create sets of installed package names (normalized)
    conda_installed = {normalize_pkg_name(pkg['name']) for pkg in installed_conda_packages}
    pip_installed = {normalize_pkg_name(pkg['name']) for pkg in installed_pip_packages}

    # Parse the environment.yml file
    try:
        with open(environment_file, 'r') as f:
            env_data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing {environment_file}: {e}")

    missing_packages = []

    # Process dependencies
    if 'dependencies' in env_data:
        for dep in env_data['dependencies']:
            if isinstance(dep, dict) and 'pip' in dep:
                # Process pip dependencies
                for pip_pkg in dep['pip']:
                    if isinstance(pip_pkg, str):
                        pkg_name = normalize_pkg_name(pip_pkg.split('==')[0].split('>=')[0].split('<=')[0].strip())
                        if pkg_name not in pip_installed and pkg_name not in conda_installed:
                            missing_packages.append(pkg_name)
            elif isinstance(dep, str):
                # Process conda dependencies
                pkg_name = normalize_pkg_name(dep.split('=')[0].strip())
                if pkg_name not in conda_installed and pkg_name not in pip_installed:
                    missing_packages.append(pkg_name)

    # Return result
    if missing_packages:
        return {
            "installed_status": False,
            "missing_packages": missing_packages
        }
    
    return {
        "installed_status": True
    }

def check_renv_lock(renv_file: str):
    """Check if R dependencies in renv.lock are installed."""
    
    # Check if the renv.lock file exists
    if not os.path.exists(renv_file):
        raise FileNotFoundError(f"{renv_file} not found.")
    
    # Get installed R packages
    try:
        installed_packages = subprocess.run(
            ["/usr/bin/Rscript", "-e", "installed.packages()[,1]"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        ).stdout.strip().split()
    except subprocess.CalledProcessError as e:
        print(f"Error getting installed R packages: {e}")
        return {
            "installed_status": False,
            "message": "Failed to get installed R packages"
        }

    # Read and parse renv.lock file
    try:
        with open(renv_file, 'r') as f:
            renv_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing renv.lock file: {e}")
        return {
            "installed_status": False,
            "message": "Invalid renv.lock file format"
        }

    # Get required packages from renv.lock
    required_packages = []
    if 'Packages' in renv_data:
        required_packages = list(renv_data['Packages'].keys())

    # Add essential packages that should be checked
    essential_packages = ['renv', 'optparse', 'gtable', 'scales', 'rlang']
    required_packages.extend(essential_packages)

    print("required_packages:", required_packages)

    # Convert installed packages to lowercase for case-insensitive comparison
    installed_packages = [pkg.lower() for pkg in installed_packages]
    
    # Check which packages are missing
    missing_packages = []
    for pkg in required_packages:
        if pkg.lower() not in installed_packages:
            missing_packages.append(pkg)

    if missing_packages:
        return {
            "installed_status": False,
            "missing_packages": missing_packages
        }
    
    return {
        "installed_status": True
    }

def verify_dependencies(dependency_file_name: str):
    """
    의존성 파일을 기반으로 시스템에 누락된 패키지를 확인합니다.

    Parameters:
        dependency_file_name (str): 의존성 파일 이름 (requirements.txt, environment.yml, renv.lock)

    Returns:
        dict: 누락된 패키지 목록과 설치 상태
    """
    # Ensure dependency_file_name is a valid string
    if not isinstance(dependency_file_name, str):
        raise ValueError(f"Invalid file name: {dependency_file_name}")

    # Check for requirements.txt
    if "requirements.txt" in dependency_file_name:
        return check_requirements_txt(dependency_file_name)
    # Check for environment.yml or environment.yaml
    elif "environment.yml" in dependency_file_name or "environment.yaml" in dependency_file_name:
        return check_environment_yml(dependency_file_name)
    # Check for renv.lock
    elif "renv.lock" in dependency_file_name:
        return check_renv_lock(dependency_file_name)
    else:
        print(f"Unsupported dependency file: {dependency_file_name}")
        return {
            "installed_status": False,
            "message": f"Unsupported file: {dependency_file_name}"
        }

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
    import tempfile
    
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

def setup_plugin_environments_via_docker(plugin_name: str):
    """
    플러그인의 Python과 R 가상환경을 설치합니다.
    """
    try:
        client = docker.DockerClient(base_url="unix://var/run/docker.sock")
        celery_container = next(
            (c for c in client.containers.list() if 'celery' in c.name),
            None
        )
        if not celery_container:
            raise RuntimeError("Celery container not found")

        plugin_path = f"/app/plugin/{plugin_name}"
        dep_path = f"{plugin_path}/dependency"
        env_path = f"{plugin_path}/env"
        py_env = f"{env_path}/py_env"
        r_env = f"{env_path}/r_env"

        def run(cmd, error_msg, workdir=None):
            result = celery_container.exec_run(cmd, workdir=workdir, user="root")
            if result.exit_code != 0:
                raise RuntimeError(f"{error_msg}: {result.output.decode()}")
            return result

        run(f"mkdir -p {env_path} {r_env}", "Failed to create environment directories")
        run(f"chmod 777 {env_path} {r_env}", "Failed to set permissions")

        plugin_path, _ = get_plugin_path(plugin_name)
        plugin_dep_path = os.path.join(plugin_path, "dependency")
        
        if os.path.exists(os.path.join(plugin_dep_path, "environment.yml")):
            run(f"micromamba create -y -p {py_env} -f {dep_path}/environment.yml", "Failed to create Python env")
        elif os.path.exists(os.path.join(plugin_dep_path, "requirements.txt")):
            run(f"micromamba create -y -p {py_env} python=3.10", "Failed to create base Python env")
            run(f"micromamba run -p {py_env} pip install -r {dep_path}/requirements.txt", "Failed to install pip packages")

        if os.path.exists(os.path.join(plugin_dep_path, "renv.lock")):
            run(f"cp {plugin_dep_path}/renv.lock {r_env}/", "Failed to copy renv.lock")
            run(f"cp {plugin_dep_path}/*.tar.gz {r_env}/ || true", "Failed to copy tar.gz files")

            r_commands = [
                "R -e 'if (!requireNamespace(\"renv\", quietly=TRUE)) install.packages(\"renv\", repos=\"https://cloud.r-project.org\")'",
                "R -e 'renv::init(bare=TRUE)'",
                "R -e 'tryCatch({ renv::restore(lockfile=\"renv.lock\", prompt=FALSE) }, error=function(e) { message(\"Warning: partial restore\") })'",
                "R -e 'writeLines(c(\"Sys.setenv(RENV_PATHS_LIBRARY=\\\"renv_library\\\")\", \"source(\\\"renv/activate.R\\\")\"), \".Rprofile\")'"
            ]
            for cmd in r_commands:
                run(cmd, "Failed to setup R environment", workdir=r_env)

        run(f"chmod -R 777 {env_path}", "Failed to finalize permissions")

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to setup plugin environment: {str(e)}")

def extract_python_version_from_environment_yml(env_path: str) -> str:
    """environment.yml 파일에서 Python 버전을 추출합니다."""
    if not os.path.isfile(env_path):
        return None
    with open(env_path, 'r') as f:
        data = yaml.safe_load(f)
    dependencies = data.get('dependencies', [])
    for dep in dependencies:
        if isinstance(dep, str) and dep.startswith("python="):
            return dep.split("=")[1]
    return None

def extract_r_version_from_renv_lock(renv_path: str) -> str:
    """renv.lock 파일에서 R 버전을 추출합니다."""
    if not os.path.isfile(renv_path):
        return None
    with open(renv_path, 'r') as f:
        data = json.load(f)
    r_info = data.get("R", {})
    if isinstance(r_info, dict):
        return r_info.get("Version", None)
    return None

def build_plugin_docker_image(plugin_path: str, plugin_name: str) -> dict:
    """
    플러그인의 Docker 이미지를 빌드하고 로그를 저장합니다.

    Parameters:
        plugin_path (str): 플러그인 폴더 경로
        plugin_name (str): 플러그인 이름

    Returns:
        dict: {
            'success': bool,
            'message': str,
            'log_file': str,
            'image_tag': str
        }
    """
    try:
        client = docker.from_env()
        client.ping()  # Docker 데몬 연결 확인
    except Exception as e:
        error_msg = f"Failed to connect to Docker daemon: {str(e)}"
        return {
            'success': False,
            'message': error_msg,
            'log_file': None,
            'image_tag': None
        }

    # 로그 디렉토리 및 파일 설정
    ensure_build_logs_dir()
    log_file = os.path.join(BUILD_LOGS_DIR, f"{plugin_name.lower()}.log")
    image_tag = f"plugin-{plugin_name.lower()}"

    # 플랫폼 자동 감지
    target_platform = get_optimal_platform(plugin_path, client)
    build_platform = get_target_platform(client)
    use_gpu = check_gpu_requirement(plugin_path)

    try:
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(f"Building Docker image for plugin: {plugin_name}\n")
            f.write(f"Start time: {datetime.now().isoformat()}\n")
            f.write(f"Docker version: {client.version().get('Version', 'unknown')}\n")
            f.write(f"Host platform: {platform.machine()} ({platform.system()})\n")
            f.write(f"Target platform: {target_platform}\n")
            f.write(f"Build platform: {build_platform}\n")
            f.write(f"GPU mode: {use_gpu}\n\n")
            f.flush()

        # 기존 이미지 제거 (있다면)
        try:
            client.images.remove(image=image_tag, force=True)
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Removed existing image: {image_tag}\n\n")
        except docker.errors.ImageNotFound:
            pass

        # Dockerfile 존재 확인
        dockerfile_path = os.path.join(plugin_path, "Dockerfile")
        if not os.path.exists(dockerfile_path):
            error_msg = f"Dockerfile not found in {plugin_path}"
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Error: {error_msg}\n")
            return {
                'success': False,
                'message': error_msg,
                'log_file': log_file,
                'image_tag': image_tag
            }

        # Dockerfile 내용 기록
        with open(log_file, "a", encoding="utf-8") as f:
            f.write("Dockerfile content:\n")
            with open(dockerfile_path, "r", encoding="utf-8") as df:
                f.write(df.read())
            f.write("\n\nStarting build...\n")
            f.flush()

        image = None
        build_output = []

        # 이미지 빌드
        try:
            # client.images.build()는 BuildKit을 올바르게 지원
            image, build_logs = client.images.build(
                path=plugin_path,
                tag=image_tag,
                rm=True,
                forcerm=True,
                buildargs={
                    'BUILDKIT_PROGRESS': 'plain',
                    'DOCKER_BUILDKIT': '1',  # BuildKit 명시적 활성화
                    'TARGETPLATFORM': target_platform,  # 자동 감지된 타겟 플랫폼
                    'BUILDPLATFORM': build_platform     # 자동 감지된 빌드 플랫폼
                },
                platform=None,  # 자동 플랫폼 감지
                pull=False,
                nocache=False
            )
            
            # 빌드 로그 처리
            for log_line in build_logs:
                with open(log_file, "a", encoding="utf-8") as f:
                    if 'stream' in log_line:
                        f.write(log_line['stream'])
                        build_output.append(log_line['stream'])
                    elif 'error' in log_line:
                        f.write(f"Build error: {log_line['error']}\n")
                        raise Exception(log_line['error'])
                    f.flush()
            
        except Exception as e:
            # Fallback: legacy API 사용 (BuildKit 없이)
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Modern API failed, trying legacy API: {str(e)}\n")
            
            for chunk in client.api.build(
                path=plugin_path,
                tag=image_tag,
                rm=True,
                forcerm=True,
                decode=True,
                buildargs={
                    'BUILDKIT_PROGRESS': 'plain',
                    'DOCKER_BUILDKIT': '0',  # Legacy API에서는 BuildKit 비활성화
                    'TARGETPLATFORM': target_platform,  # 자동 감지된 타겟 플랫폼
                    'BUILDPLATFORM': build_platform     # 자동 감지된 빌드 플랫폼
                }
            ):
                with open(log_file, "a", encoding="utf-8") as f:
                    if 'stream' in chunk:
                        f.write(chunk['stream'])
                        build_output.append(chunk['stream'])
                    elif 'error' in chunk:
                        f.write(f"Build error: {chunk['error']}\n")
                        raise Exception(chunk['error'])
                    elif 'aux' in chunk and 'ID' in chunk['aux']:
                        image = client.images.get(chunk['aux']['ID'])
                    f.flush()

        if not image:
            image = client.images.get(image_tag)

        # 빌드 성공 기록
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"\nBuild completed successfully at {datetime.now().isoformat()}\n")
            if image:
                f.write(f"Image ID: {image.id}\n")
                f.write(f"Size: {image.attrs['Size'] / 1024 / 1024:.2f} MB\n")
                f.write(f"Created: {image.attrs['Created']}\n")
            f.flush()

        return {
            'success': True,
            'message': f"Successfully built image: {image_tag}",
            'log_file': log_file,
            'image_tag': image_tag
        }

    except Exception as e:
        error_msg = f"Unexpected error during Docker build: {str(e)}"
        try:
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"\nCritical Error: {error_msg}\n")
        except:
            pass
        return {
            'success': False,
            'message': error_msg,
            'log_file': log_file,
            'image_tag': image_tag
        }

def generate_base_image_section(use_gpu=True):
    """Base image section generation"""
    lines = []
    if use_gpu:
        lines.append("FROM --platform=linux/amd64 nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04")
    else:
        lines.append("ARG TARGETPLATFORM=linux/amd64")
        lines.append("ARG BUILDPLATFORM=linux/amd64")
        lines.append("")
        lines.append("FROM --platform=$TARGETPLATFORM debian:bullseye-slim")
    return lines


def generate_env_setup_section():
    """Environment variables setup section"""
    return [
        "",
        "ENV DEBIAN_FRONTEND=noninteractive",
        "ENV TZ=Asia/Seoul"
    ]


def generate_system_packages_section(use_gpu=True, has_r=False):
    """Enhanced system packages installation section based on Scribe requirements"""
    lines = [
        "",
        "RUN apt-get update && apt-get install -y \\"
    ]
    
    # Common packages
    common_packages = [
        "    build-essential gcc g++ gfortran make \\",
        "    libssl-dev libcurl4-openssl-dev libxml2-dev libxslt1-dev \\",
        "    libjpeg-dev libpng-dev libfreetype6-dev libtiff5-dev \\",
        "    libx11-dev libxt-dev" + (" xorg-dev \\" if has_r else " \\"),
        "    libglu1-mesa-dev \\",
        "    libharfbuzz-dev libfribidi-dev \\",
        "    libglpk-dev \\",
        "    curl wget unzip git \\",
        "    python3 python3-pip python3-venv \\",
        "    software-properties-common \\",
        "    gnupg ca-certificates \\",
        "    zlib1g-dev" + (" liblzma-dev \\" if has_r else "")
    ]
    
    lines.extend(common_packages)
    
    # Comprehensive R ecosystem packages (optimized for Scribe)
    if has_r:
        r_packages = [
            "    rsync \\",
            "    libv8-dev libnode-dev \\",
            "    libudunits2-dev libgdal-dev \\",
            "    libcairo2-dev \\",
            "    libfontconfig1-dev \\", 
            "    libgeos-dev \\",
            "    libproj-dev \\",
            "    libmagick++-dev \\",
            "    libpoppler-cpp-dev \\",
            "    librsvg2-dev \\"
        ]
        lines.extend(r_packages)
    
    # 마지막 패키지 라인이 백슬래시로 끝나지 않으면 백슬래시 추가
    if not lines[-1].endswith(" \\"):
        lines[-1] = lines[-1] + " \\"
    
    lines.append("    && apt-get clean && rm -rf /var/lib/apt/lists/*")
    return lines


def generate_dependency_copy_section():
    """Dependency folder copy section"""
    return [
        "",
        "COPY dependency/ /workspace/dependency/"
    ]


def generate_micromamba_install_section(use_gpu=False):
    """Micromamba installation section - dynamic platform-specific installation"""
    if use_gpu:
        # GPU 버전은 linux/amd64 고정
        return [
            "",
            "RUN mkdir -p /usr/local/bin && \\",
            "    cd /tmp && \\",
            "    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj && \\",
            "    cp bin/micromamba /usr/local/bin/micromamba && \\",
            "    chmod +x /usr/local/bin/micromamba && \\",
            "    rm -rf /tmp/bin /tmp/info && \\",
            "    /usr/local/bin/micromamba --version"
        ]
    else:
        # CPU 버전은 플랫폼별 동적 설치 (기본값 포함)
        return [
            "",
            "RUN PLATFORM=${TARGETPLATFORM:-linux/amd64} && \\",
            "    case $PLATFORM in \\",
            "        \"linux/amd64\")  MICROMAMBA_ARCH=linux-64  ;; \\",
            "        \"linux/arm64\")  MICROMAMBA_ARCH=linux-aarch64  ;; \\",
            "        \"linux/arm/v7\") MICROMAMBA_ARCH=linux-armv7l  ;; \\",
            "        *) echo \"Unsupported platform: $PLATFORM, defaulting to linux-64\" && MICROMAMBA_ARCH=linux-64 ;; \\",
            "    esac && \\",
            "    mkdir -p /usr/local/bin && \\",
            "    curl -Ls \"https://micro.mamba.pm/api/micromamba/${MICROMAMBA_ARCH}/latest\" | tar -xvj -C /tmp && \\",
            "    cp /tmp/bin/micromamba /usr/local/bin/micromamba && \\",
            "    chmod +x /usr/local/bin/micromamba && \\",
            "    rm -rf /tmp/bin /tmp/info && \\",
            "    /usr/local/bin/micromamba --version"
        ]


def generate_env_variables_section(has_r=False):
    """Environment variables configuration section"""
    lines = [
        "",
        "ENV MAMBA_ROOT_PREFIX=/opt/micromamba",
        "ENV MAMBA_EXE=/usr/local/bin/micromamba",
        "ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$MAMBA_ROOT_PREFIX/bin:$PATH"
    ]
    
    if has_r:
        lines.append("ENV R_HOME=/opt/micromamba/envs/plugin_env/lib/R")
    
    return lines


def generate_micromamba_setup_section():
    """Micromamba environment directory creation section"""
    return [
        "",
        "RUN mkdir -p $MAMBA_ROOT_PREFIX && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/envs && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/pkgs && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/etc/profile.d"
    ]


def generate_python_env_section(python_version="3.10"):
    """Python environment creation section"""
    return [
        "",
        f"RUN /usr/local/bin/micromamba create -y -n plugin_env python={python_version} -c conda-forge --root-prefix $MAMBA_ROOT_PREFIX"
    ]


def generate_snakemake_install_section():
    """Snakemake and essential packages installation section"""
    return [
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    pip install --no-cache-dir \\",
        "    'snakemake==7.14.0' \\",
        "    'pulp==2.7.0' \\",
        "    'tabulate==0.8.10'"
    ]


def generate_python_packages_section(has_requirements=False):
    """Python packages installation section"""
    if not has_requirements:
        return []
    
    return [
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    pip install --no-cache-dir -r /workspace/dependency/requirements.txt || true"
    ]


def generate_r_installation_section(r_version="4.4.2", has_renv=False, use_gpu=False):
    """Enhanced R installation section - efficient package management with rsync"""
    lines = []
    
    # Comprehensive Micromamba R packages installation
    lines.extend([
        "",
        f"RUN /usr/local/bin/micromamba install -y -n plugin_env \\",
        f"    r-base={r_version} \\",
        "    r-essentials \\",
        "    r-renv \\",
        "    r-jsonlite \\",
        "    zlib \\",
        "    xz \\",
        "    libpng \\",
        "    libjpeg-turbo \\",
        "    libcurl \\",
        "    libxml2 \\",
        "    libxslt \\",
        "    cairo \\",
        "    pango \\",
        "    harfbuzz \\",
        "    fribidi \\",
        "    freetype \\",
        "    fontconfig \\",
        "    libtiff \\",
        "    pkg-config \\",
        "    -c conda-forge \\",
        "    -r $MAMBA_ROOT_PREFIX && \\",
        "    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/R /usr/local/bin/R && \\",
        "    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/Rscript /usr/local/bin/Rscript"
    ])
    
    # R configuration
    lines.extend([
        "",
        "RUN echo 'options(repos = c(CRAN = \"https://cloud.r-project.org\"), download.file.method = \"libcurl\")' > /root/.Rprofile && \\",
        "    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    Rscript -e \"install.packages(c('renv', 'BiocManager'), dependencies = TRUE)\""
    ])
    
    # Enhanced renv packages installation
    if has_renv:
        lines.extend([
            "",
            "RUN if [ -f \"/workspace/dependency/renv.lock\" ]; then \\",
            "    echo 'Installing R packages from renv.lock...' && \\",
            "    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
            "    Rscript -e \"\\",
            "        # Set micromamba environment as primary library path \\",
            "        target_lib <- '/opt/micromamba/envs/plugin_env/lib/R/library'; \\",
            "        .libPaths(c(target_lib, .libPaths())); \\",
            "        cat('Library paths configured:', paste(.libPaths(), collapse=', '), '\\\\n'); \\",
            "        \\",
            "        # Change to dependency directory for renv operations \\",
            "        setwd('/workspace/dependency'); \\",
            "        \\",
            "        # Load renv and configure \\",
            "        library(renv); \\",
            "        options(renv.config.cache.enabled = FALSE); \\",
            "        options(renv.config.auto.snapshot = FALSE); \\",
            "        \\",
            "        # Initialize renv project \\",
            "        cat('Initializing renv project...\\\\n'); \\",
            "        renv::init(bare = TRUE, restart = FALSE); \\",
            "        \\",
            "        # Restore packages to renv library \\",
            "        cat('Restoring packages from renv.lock...\\\\n'); \\",
            "        tryCatch({ \\",
            "            renv::restore(lockfile = '/workspace/dependency/renv.lock', confirm = FALSE); \\",
            "            cat('✅ R packages successfully installed to renv library\\\\n'); \\",
            "        }, error = function(e) { \\",
            "            cat('⚠️ Some packages may have failed, but continuing...\\\\n'); \\",
            "            cat('Error details:', conditionMessage(e), '\\\\n'); \\",
            "        }); \\",
            "    \" && \\",
            "    echo 'renv installation phase completed' ; \\",
            "fi",
            "",
            "RUN if [ -d \"/workspace/dependency/renv/library\" ]; then \\",
            "    echo '=== Stage 1: Starting complete rsync copy ===' && \\",
            "    RENV_LIB_PATH=\"/workspace/dependency/renv/library/linux-debian-bullseye/R-4.4/x86_64-conda-linux-gnu\" && \\",
            "    TARGET_LIB_PATH=\"/opt/micromamba/envs/plugin_env/lib/R/library\" && \\",
            "    if [ -d \"$RENV_LIB_PATH\" ]; then \\",
            "        echo \"Copy source: $RENV_LIB_PATH\" && \\",
            "        echo \"Copy target: $TARGET_LIB_PATH\" && \\",
            "        echo \"Packages to copy:\" && \\",
            "        ls -1 \"$RENV_LIB_PATH\" | wc -l && echo \" packages\" && \\",
            "        echo \"Starting rsync copy...\" && \\",
            "        rsync -av --progress \"$RENV_LIB_PATH/\" \"$TARGET_LIB_PATH/\" && \\",
            "        echo \"✅ rsync copy completed\" && \\",
            "        echo \"Post-copy package verification:\" && \\",
            "        ls -1 \"$TARGET_LIB_PATH\" | wc -l && echo \" packages\" ; \\",
            "    else \\",
            "        echo \"⚠️ Cannot find renv library path: $RENV_LIB_PATH\" ; \\",
            "    fi ; \\",
            "else \\",
            "    echo \"⚠️ renv library directory does not exist\" ;"
            "fi",
            "",
            "RUN echo '=== Stage 2: Cleanup temporary files for space optimization ===' && \\",
            "    echo '=== Stage 2: Starting temporary file cleanup ===' && \\",
            "    # Remove renv temporary library directory \\",
            "    if [ -d \"/workspace/dependency/renv/library\" ]; then \\",
            "        echo 'Removing renv temporary library...' && \\",
            "        RENV_LIB_SIZE=$(du -sh /workspace/dependency/renv/library 2>/dev/null | cut -f1 || echo \"unknown\") && \\",
            "        echo \"renv library size to remove: $RENV_LIB_SIZE\" && \\",
            "        rm -rf /workspace/dependency/renv/library && \\",
            "        echo '✅ renv temporary library removal completed' ; \\",
            "    fi && \\",
            "    # Remove renv cache \\",
            "    if [ -d \"/root/.cache/R/renv\" ]; then \\",
            "        echo 'Removing renv cache...' && \\",
            "        CACHE_SIZE=$(du -sh /root/.cache/R/renv 2>/dev/null | cut -f1 || echo \"unknown\") && \\",
            "        echo \"Cache size to remove: $CACHE_SIZE\" && \\",
            "        rm -rf /root/.cache/R/renv && \\",
            "        echo '✅ renv cache removal completed' ; \\",
            "    fi && \\",
            "    # Remove rsync package \\",
            "    echo 'Removing rsync package...' && \\",
            "    apt-get remove -y rsync && \\",
            "    apt-get autoremove -y && \\",
            "    echo '✅ rsync package removal completed' && \\",
            "    echo '=== Stage 2: Temporary file cleanup completed ==='"
        ])
    
    lines.extend([
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    Rscript -e \"if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', dependencies = TRUE)\"",
        "",
        "ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE"
    ])
    
    return lines

def generate_workspace_setup_section():
    """Workspace setup section"""
    return [
        "",
        "RUN mkdir -p /workspace/logs && \\",
        "    chmod 777 /workspace"
    ]


def generate_copy_files_section(plugin_path):
    """Snakefile and scripts copy section"""
    lines = [
        "",
        "COPY Snakefile /workspace/Snakefile",
        "",
        "COPY scripts/ /scripts/",
        "",
        "WORKDIR /workspace"
    ]
    
    return lines


def generate_entrypoint_section(has_r=False):
    """Entrypoint script generation section - enhanced version"""
    lines = [
        "",
        "RUN echo '#!/bin/bash' > /entrypoint.sh && \\",
        "    echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> /entrypoint.sh && \\",
        "    echo 'export MAMBA_EXE=/usr/local/bin/micromamba' >> /entrypoint.sh && \\",
        "    echo 'export PATH=$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$PATH' >> /entrypoint.sh && \\"
    ]
    
    if has_r:
        lines.extend([
            "    echo 'export R_HOME=/opt/micromamba/envs/plugin_env/lib/R' >> /entrypoint.sh && \\",
            "    echo 'export R_LIBS_USER=/opt/micromamba/envs/plugin_env/lib/R/library' >> /entrypoint.sh && \\",
            "    echo 'export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /entrypoint.sh && \\"
        ])
    
    lines.extend([
        "    echo 'eval \"$($MAMBA_EXE shell activate -s bash -p $MAMBA_ROOT_PREFIX plugin_env)\" 2>/dev/null || true' >> /entrypoint.sh && \\",
        "    echo 'cd /workspace' >> /entrypoint.sh && \\",
        "    echo 'exec \"$@\"' >> /entrypoint.sh && \\",
        "    chmod +x /entrypoint.sh",
        "",
        "ENTRYPOINT [\"/entrypoint.sh\"]"
    ])
    
    return lines


def generate_healthcheck_section():
    """Health check section"""
    return [
        "",
        "HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \\",
        "    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1"
    ]


def generate_cmd_section():
    """Default command section"""
    return [
        "",
        "CMD [\"/bin/bash\"]"
    ]


def generate_plugin_dockerfile(plugin_path: str, output_path: str, use_gpu: bool = True):
    """
    Analyzes plugin folder and generates multi-platform Dockerfile.
    
    Parameters:
        plugin_path (str): Plugin folder path
        output_path (str): Output Dockerfile path
        use_gpu (bool): GPU usage flag (default: True)
    """
    dependency_path = os.path.join(plugin_path, "dependency")
    
    # Check dependency files
    has_requirements = os.path.isfile(os.path.join(dependency_path, "requirements.txt"))
    has_environment = os.path.isfile(os.path.join(dependency_path, "environment.yml"))
    has_renv = os.path.isfile(os.path.join(dependency_path, "renv.lock"))
    
    # Check Python and R scripts
    has_python = False
    has_r = False
    scripts_path = os.path.join(plugin_path, "scripts")
    if os.path.exists(scripts_path):
        for file in os.listdir(scripts_path):
            if file.endswith(".py"):
                has_python = True
            if file.endswith(".R"):
                has_r = True
    
    # Extract versions
    python_version = "3.10"  # Default
    r_version = "4.4.2"  # Default
    
    if has_environment:
        python_version = extract_python_version_from_environment_yml(
            os.path.join(dependency_path, "environment.yml")
        )
    
    if has_renv:
        r_version = extract_r_version_from_renv_lock(
            os.path.join(dependency_path, "renv.lock")
        )
    
    # Generate Dockerfile
    dockerfile_lines = []
    
    # 1. Base image (multi-platform support)
    dockerfile_lines.extend(generate_base_image_section(use_gpu))
    
    # 2. Environment variables setup
    dockerfile_lines.extend(generate_env_setup_section())
    
    # 3. System packages installation
    dockerfile_lines.extend(generate_system_packages_section(use_gpu, has_r or has_renv))
    
    # 4. Dependency folder copy
    dockerfile_lines.extend(generate_dependency_copy_section())
    
    # 5. Micromamba installation (platform-specific)
    dockerfile_lines.extend(generate_micromamba_install_section(use_gpu))
    
    # 6. Environment variables configuration
    dockerfile_lines.extend(generate_env_variables_section(has_r or has_renv))
    
    # 7. Micromamba environment directory creation
    dockerfile_lines.extend(generate_micromamba_setup_section())
    
    # 8. Python environment creation
    dockerfile_lines.extend(generate_python_env_section(python_version))
    
    # 9. Snakemake and essential packages installation
    dockerfile_lines.extend(generate_snakemake_install_section())
    
    # 10. Python packages installation
    if has_requirements:
        dockerfile_lines.extend(generate_python_packages_section(has_requirements))
    
    # 11. R installation and packages (platform consideration)
    if has_r or has_renv:
        dockerfile_lines.extend(generate_r_installation_section(r_version, has_renv, use_gpu))
    
    # 12. Workspace setup
    dockerfile_lines.extend(generate_workspace_setup_section())
    
    # 13. Snakefile and scripts copy
    dockerfile_lines.extend(generate_copy_files_section(plugin_path))
    
    # 14. Entrypoint configuration
    dockerfile_lines.extend(generate_entrypoint_section(has_r or has_renv))
    
    # 15. Health check configuration
    dockerfile_lines.extend(generate_healthcheck_section())
    
    # 16. Default command configuration
    dockerfile_lines.extend(generate_cmd_section())
    
    # Save Dockerfile
    with open(output_path, "w") as f:
        f.write("\n".join(dockerfile_lines))
    
    # Output result logs
    print(f"[✓] Multi-platform Dockerfile generated for {plugin_path}")
    print(f"    - GPU: {use_gpu}")
    print(f"    - Platform support: {'linux/amd64' if use_gpu else 'linux/amd64, linux/arm64, linux/arm/v7'}")
    print(f"    - Python version: {python_version}")
    if has_r or has_renv:
        print(f"    - R version: {r_version}")
    print(f"    - has_requirements: {has_requirements}")
    print(f"    - has_renv: {has_renv}")

def check_plugin_docker_image(plugin_name: str) -> bool:
    """
    플러그인의 Docker 이미지가 존재하는지 확인합니다.

    Parameters:
        plugin_name (str): 플러그인 이름

    Returns:
        bool: Docker 이미지 존재 여부
    """
    try:
        client = docker.from_env()
        image_tag = f"plugin-{plugin_name.lower()}"
        
        # 이미지 존재 여부 확인
        try:
            client.images.get(image_tag)
            return True
        except docker.errors.ImageNotFound:
            return False
        except Exception as e:
            print(f"Error checking Docker image: {str(e)}")
            return False
            
    except Exception as e:
        print(f"Failed to connect to Docker daemon: {str(e)}")
        return False

def generate_plugin_image_uri(plugin_name: str, source: str, version: Optional[str] = None) -> str:
    """
    Generate plugin image URI based on plugin type.
    
    Args:
        plugin_name: Name of the plugin
        source: Plugin source ("official" or "local")  
        version: Version for official plugins (optional)
        
    Returns:
        Plugin image URI string
    """
    if source == "official":
        # For official plugins, use GitHub Container Registry
        version = version or "latest"
        return f"ghcr.io/cxinsys/cellcraft-{plugin_name.lower()}:{version}"
    else:
        # For local plugins, use timestamp-based identifier
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return f"local-build:{plugin_name.lower()}-{timestamp}"
