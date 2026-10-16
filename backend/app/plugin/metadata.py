"""Plugin metadata & descriptor helpers (split from ``plugin/utils.py`` in PR-9).

Two responsibilities that both describe a plugin rather than build/run it:

- Dependency manifest verification (requirements.txt / environment.yml /
  renv.lock) via ``verify_dependencies`` and its per-format checkers.
- Drawflow workflow-graph template generation for a plugin's default node
  layout (``generate_plugin_drawflow_template`` /
  ``generate_merged_plugin_drawflow``).

The drawflow generators are placed here (rather than a dedicated module) because
they describe the plugin's default workflow representation and the 800-line
per-file limit left ``metadata.py`` as the only module with capacity for them.
Function signatures and behavior are unchanged from the original ``utils.py``.
"""
import os
import json
import subprocess

import yaml
from typing import List, Dict, Any


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
    pos_x_resultfile = pos_x_algorithm + pos_x_increment
    pos_x_visualization = pos_x_resultfile + pos_x_increment

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
    algorithm_node_id = node_id
    node_id += 1

    # ResultFiles 노드 생성
    resultfiles_node_id = node_id
    resultfiles_node = {
        "id": resultfiles_node_id,
        "name": "ResultFiles",
        "data": {"title": "ResultFiles"},
        "class": "ResultFiles",
        "html": "ResultFiles",
        "typenode": "vue",
        "inputs": {
            "input_1": {
                "connections": [{"node": str(algorithm_node_id), "input": "output_1"}]
            }
        },
        "outputs": {
            "output_1": {
                "connections": []
            }
        },
        "pos_x": pos_x_resultfile,
        "pos_y": pos_y_start  # Same Y as Algorithm node
    }
    merged_drawflow["drawflow"]["Home"]["data"][str(node_id)] = resultfiles_node
    node_id += 1

    # Visualization 노드 생성
    visualization_node = {
        "id": node_id,
        "name": "Visualization",
        "data": {},  # No title
        "class": "Visualization",
        "html": "Visualization",
        "typenode": "vue",
        "inputs": {
            "input_1": {
                "connections": [{"node": str(resultfiles_node_id), "input": "output_1"}]
            }
        },
        "outputs": {},  # No outputs
        "pos_x": pos_x_visualization,
        "pos_y": pos_y_start  # Same Y as Algorithm node
    }
    merged_drawflow["drawflow"]["Home"]["data"][str(node_id)] = visualization_node

    # Update Algorithm output connection
    algorithm_node["outputs"]["output_1"]["connections"].append({
        "node": str(resultfiles_node_id),
        "output": "input_1"
    })

    # Update ResultFiles output connection
    resultfiles_node["outputs"]["output_1"]["connections"].append({
        "node": str(visualization_node['id']),
        "input": "input_1"
    })

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

    # ResultFiles 노드 생성
    resultfiles_node_id = node_id
    resultfiles_node = {
        "id": resultfiles_node_id,
        "name": "ResultFiles",
        "data": {"title": "ResultFiles"},
        "class": "ResultFiles",
        "html": "ResultFiles",
        "typenode": "vue",
        "inputs": {
            "input_1": {
                "connections": [{"node": str(algorithm_node_id), "input": "output_1"}]
            }
        },
        "outputs": {
            "output_1": {
                "connections": []
            }
        },
        "pos_x": pos_x_resultfile,
        "pos_y": pos_y_algorithm  # Same Y as Algorithm node
    }
    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = resultfiles_node
    node_id += 1

    # Visualization 노드 생성
    visualization_node = {
        "id": node_id,
        "name": "Visualization",
        "data": {},  # No title
        "class": "Visualization",
        "html": "Visualization",
        "typenode": "vue",
        "inputs": {
            "input_1": {
                "connections": [{"node": str(resultfiles_node_id), "input": "output_1"}]
            }
        },
        "outputs": {},  # No outputs
        "pos_x": pos_x_visualization,
        "pos_y": pos_y_algorithm  # Same Y as Algorithm node
    }
    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = visualization_node
    node_id += 1

    # Update Algorithm output connection
    algorithm_node["outputs"]["output_1"]["connections"].append({
        "node": str(resultfiles_node_id),
        "output": "input_1"
    })

    # Update ResultFiles output connection
    resultfiles_node["outputs"]["output_1"]["connections"].append({
        "node": str(visualization_node['id']),
        "input": "input_1"
    })

    # 생성된 new_drawflow 데이터를 출력합니다.
    print(json.dumps(new_drawflow, indent=2))
    return new_drawflow
