import json
import os
import re
import subprocess
import yaml
import shutil
from typing import Dict, Any
import numpy as np
from fastapi import HTTPException

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
            if not output_val['connections']:
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

def generate_snakemake_code(rules_data, output_folder_path):
    snakemake_path = os.path.join(output_folder_path, "Snakefile")
    visualization_snakemake_path = os.path.join(output_folder_path, "visualization_Snakefile")

    snakemake_code = "import os\n\n"  # Add import statement at the top
    visualization_snakemake_code = "import os\n\n"
    input_output_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
    logs_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/logs"
    unique_input_path = "user/{user_name}/data"

    # Find all outputs across all rules to determine unique inputs not present in outputs
    all_outputs = {out for rule in rules_data.values() for out in rule['output']}
    all_inputs = {inp for rule in rules_data.values() for inp in rule['input']}
    unique_inputs = [inp for inp in all_inputs if inp not in all_outputs]
    scripts_log_word = "> {log.stdout} 2> {log.stderr}"

    # Iterate through each rule in the dictionary
    for rule_id, rule in rules_data.items():
        # Determine which code block to append to: snakemake_code or visualization_snakemake_code
        if rule.get('isVisualization', False):
            target_code = visualization_snakemake_code
            print(f"Adding to visualization_snakemake_code: {rule['name']}")
        else:
            target_code = snakemake_code
            print(f"Adding to snakemake_code: {rule['name']}")
        
        # Start rule block with rule name
        target_code += f"rule {normalize_param_name(rule['name'])}:\n"

        # Input section with optional input handling
        if 'input' in rule and rule['input']:
            def get_input_param_name(inp, rule_params):
                for param in rule_params:
                    if param['type'] == 'inputFile' and param.get('defaultValue') == inp:
                        return normalize_param_name(param['name'])
                return normalize_param_name(os.path.splitext(os.path.basename(inp))[0])

            input_files = ",\n        ".join([
                f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_output_path}/{{{inp}}}\"" 
                if 'target' in get_input_param_name(inp, rule.get('parameters', [])) else
                f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{unique_input_path}/{{{inp}}}\"" 
                if inp in unique_inputs else
                f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_output_path}/{inp}\""
                for inp in rule['input']
                if not any(param['type'] == 'optionalInputFile' and param.get('defaultValue') == inp 
                          for param in rule.get('parameters', []))
            ])
            if input_files:  # Only add input section if there are non-optional input files
                target_code += f"    input:\n        {input_files}\n"

        # Output section
        if 'output' in rule and rule['output']:
            def get_output_param_name(out, rule_params):
                for param in rule_params:
                    if param['type'] == 'outputFile' and param.get('defaultValue') == out:
                        return normalize_param_name(param['name'])
                return normalize_param_name(os.path.splitext(os.path.basename(out))[0])

            if rule.get('isVisualization', False):
                additional_path = "{visualization_result_path}"
                output_files = ",\n        ".join([
                    f"{get_output_param_name(out, rule.get('parameters', []))}=\"{input_output_path}/{additional_path}{out}\""
                    for out in rule['output']
                ])
            else:
                output_files = ",\n        ".join([
                    f"{get_output_param_name(out, rule.get('parameters', []))}=\"{input_output_path}/{out}\""
                    for out in rule['output']
                ])
            target_code += f"    output:\n        {output_files}\n"

        # Params section
        if 'parameters' in rule and rule['parameters']:
            param_list = []
            for param in rule['parameters']:
                normalized_name = normalize_param_name(param['name'])
                if param['name'] == "clusters" and param['type'] == "h5adParameter":
                    param_list.append(f'clusters=lambda wc: ";".join({{{param["name"]}}})')
                elif param['type'] == 'optionalInputFile':
                    param_list.append(
                        f"{normalized_name}=lambda wildcards: \"{unique_input_path}/{{{param['defaultValue']}}}\" "
                        f"if os.path.exists(\"{unique_input_path}/{{{param['defaultValue']}}}\") else \"None\""
                    )
                elif param['type'] == 'string' and param['name'] == "ScatterPlot":
                    param_list.append(
                        'ScatterPlot=lambda wildcards: {ScatterPlot} if os.path.exists({ScatterPlot}) else "None"'
                    )
                elif param['type'] != 'inputFile' and param['type'] != 'outputFile':
                    param_list.append(f"{normalized_name}={{{param['name']}}}")

            if param_list:
                param_list_str = ",\n        ".join(param_list)
                target_code += f"    params:\n        {param_list_str}\n"

        # Log section
        target_code += f"    log:\n"
        target_code += f"        stdout=\"{logs_path}/{normalize_param_name(rule['name'])}.stdout\",\n"
        target_code += f"        stderr=\"{logs_path}/{normalize_param_name(rule['name'])}.stderr\"\n"

        # Shell section based on script type
        if 'script' in rule and rule['script']:
            script_path = "/app/plugin/{plugin_name}/scripts/" + normalize_string(rule['script'])
            if rule['script'].endswith('.py'):
                shell_command = f"/opt/conda/envs/snakemake/bin/python {script_path}"
            elif rule['script'].endswith('.R'):
                shell_command = f"/usr/bin/Rscript {script_path}"
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

            target_code += f"    shell:\n        \"{shell_command} {scripts_log_word}\"\n"

        target_code += "\n"  # Add a newline for separation between rules

        # Assign modified code back to the correct variable
        if rule.get('isVisualization', False):
            visualization_snakemake_code = target_code
        else:
            snakemake_code = target_code

    # Write the generated Snakemake code to a file
    with open(snakemake_path, 'w') as file:
        file.write(snakemake_code)

    # Write the generated visualization Snakefile
    with open(visualization_snakemake_path, 'w') as file:
        file.write(visualization_snakemake_code)

    print(f"Snakemake code has been generated and saved to {output_folder_path}.")

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
    Verify if dependencies from the given file are correctly installed.
    
    Parameters:
        dependency_file_name (str): The name of the dependency file (requirements.txt, environment.yml, renv.lock)
    
    Returns:
        dict: Installed status and missing packages if any.
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

def install_dependencies(dependency_file_name: str):
    """
    의존성 파일을 기반으로 시스템에 누락된 패키지만 설치합니다.

    Parameters:
        dependency_file_name (str): 의존성 파일 이름 (requirements.txt, environment.yml, renv.lock)
    """
    dependency_file_path = os.path.abspath(dependency_file_name)

    # 의존성 검사 수행
    check_result = verify_dependencies(dependency_file_path)
    if check_result.get("installed_status", False):
        print("모든 의존성이 이미 설치되어 있습니다.")
        return

    missing_packages = check_result.get("missing_packages", [])
    if not missing_packages:
        print("설치할 의존성이 없습니다.")
        return

    print(f"설치가 필요한 패키지들: {', '.join(missing_packages)}")

    # 실행 환경 변수 설정
    conda_env_path = "/opt/conda/envs/snakemake"
    pip_executable = f"{conda_env_path}/bin/pip"
    conda_executable = f"{conda_env_path}/bin/conda"
    rscript_executable = "/usr/bin/Rscript"

    try:
        if dependency_file_name.endswith("requirements.txt"):
            print(f"{dependency_file_path}에서 Python 의존성 설치 중...")
            subprocess.run([pip_executable, "install", "--cache-dir", "/tmp/pip_cache", "-r", dependency_file_path], check=True)

        elif dependency_file_name.endswith(("environment.yml", "environment.yaml")):
            print(f"{dependency_file_path}에서 Conda 의존성 설치 중...")
            subprocess.run([conda_executable, "env", "update", "--file", dependency_file_path, "--prune"], check=True)

        elif dependency_file_name.endswith("renv.lock"):
            print(f"{dependency_file_path}에서 R 의존성 설치 중...")
            dependency_dir = os.path.dirname(dependency_file_path)

            # renv.lock 파일 로드
            with open(dependency_file_path, "r") as f:
                renv_data = json.load(f)

            local_packages = {}
            repo_local_packages = {}

            for pkg, data in renv_data.get("Packages", {}).items():
                if data.get("Source") == "Local" and os.path.exists(os.path.join(dependency_dir, data["Path"])):
                    local_packages[pkg] = data["Path"]
                elif data.get("Source") == "Repository" and data.get("Repository") == "Local" and os.path.exists(os.path.join(dependency_dir, data["Path"])):
                    repo_local_packages[pkg] = data["Path"]

            # R 설치 스크립트 실행
            renv_commands = [
                f'setwd("{dependency_dir}");',
                'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos="https://cloud.r-project.org");',
                'options(download.file.method = "libcurl");',
                'Sys.setenv(R_HOME_USER = "/root");',  # Docker 환경 대응
                'Sys.setenv(R_LIBS_USER = "/usr/local/lib/R/site-library");',  # 패키지 저장 위치 설정
                'renv::settings$use.cache(FALSE);',
                'renv::repair();',  # 패키지 손상 복구
                'renv::restore(lockfile = "renv.lock", prompt = FALSE);'
            ]

            # 기본 패키지 (`optparse`, `gtable`, `scales`, `rlang`) 추가 설치
            base_packages = ["optparse", "gtable", "scales", "rlang"]
            renv_commands.append(
                'missing_base_packages <- setdiff(c("optparse", "gtable", "scales", "rlang"), rownames(installed.packages())); '
                'if (length(missing_base_packages) > 0) install.packages(missing_base_packages, repos="https://cloud.r-project.org");'
            )

            # Local 패키지 선행 설치 (e.g., RANNinf)
            for pkg, path in local_packages.items():
                renv_commands.append(f'renv::install(local("{path}"));')

            # Repository 기반 Local 패키지 후행 설치 (e.g., Scribe)
            for pkg, path in repo_local_packages.items():
                renv_commands.append(f'renv::install(local("{path}"));')

            # R 스크립트 실행
            subprocess.run([rscript_executable, "-e", " ".join(renv_commands)], check=True)

        else:
            print(f"지원하지 않는 파일 형식입니다: {dependency_file_name}")

    except subprocess.CalledProcessError as e:
        print(f"의존성 설치 중 오류 발생: {str(e)}")
        raise

    print("누락된 의존성 설치가 완료되었습니다.")

def create_plugin_folder(plugin_folder: str):
    """
    Create the plugin folder if it doesn't exist.
    
    Parameters:
        plugin_folder (str): Path to the plugin folder.
    """
    if not os.path.exists(plugin_folder):
        os.makedirs(plugin_folder)
        print(f"Plugin folder created at {plugin_folder}")
    else:
        print(f"Plugin folder already exists at {plugin_folder}")

def create_dependency_folder(dependency_folder: str, dependencies: dict):
    """
    Create the dependency folder and add dependency files.

    Parameters:
        dependency_folder (str): Path to the dependency folder.
        dependencies (dict): A dictionary where the keys are file names and values are file contents.
    
    Raises:
        HTTPException: If the format of dependencies is invalid.
    """
    # Check if the dependency folder exists, create it if necessary
    if not os.path.exists(dependency_folder):
        os.makedirs(dependency_folder)
        print(f"Dependency folder created at {dependency_folder}")
    else:
        # dependency 폴더 안에 있는 파일들 삭제
        for file in os.listdir(dependency_folder):
            file_path = os.path.join(dependency_folder, file)
            os.remove(file_path)
        print(f"Dependency folder already exists at {dependency_folder}")

    # Validate dependencies format and create files
    if dependencies is None:
        dependencies = {}
    elif isinstance(dependencies, dict):
        for file_name, file_content in dependencies.items():
            dep_path = os.path.join(dependency_folder, file_name)
            with open(dep_path, 'w') as f:
                f.write(file_content)
            print(f"Dependency file created: {dep_path}")
    else:
        raise HTTPException(status_code=400, detail="Invalid dependencies format")


def create_reference_folder(script_folder: str, reference_folders: dict):
    """
    Create folders and files in the reference folder structure.

    Parameters:
        script_folder (str): Base path where folders and files will be created.
        reference_folders (dict): A nested dictionary representing folder and file structures.

    Raises:
        HTTPException: If the format of reference_folders is invalid.
    """

    # Check if the script_folder exists, and clear its contents if necessary
    if os.path.exists(script_folder):
        shutil.rmtree(script_folder)  # Remove entire folder and its contents
        print(f"Script folder '{script_folder}' cleared.")

    os.makedirs(script_folder)
    print(f"Script folder '{script_folder}' created.")

    # Helper function to process folders recursively
    def process_folder(current_path: str, folder_data: dict):
        for name, content in folder_data.items():
            if name == "subFolders" and isinstance(content, list):
                # subFolders 리스트가 있을 경우, 재귀적으로 처리
                for sub_folder in content:
                    for sub_folder_name, sub_folder_data in sub_folder.items():
                        folder_path = os.path.join(current_path, sub_folder_name)
                        if not os.path.exists(folder_path):
                            os.makedirs(folder_path, exist_ok=True)
                            print(f"Folder created: {folder_path}")
                        process_folder(folder_path, sub_folder_data)  # Recur into the subfolder
            elif isinstance(content, dict):  # 일반적인 서브 폴더 처리
                folder_path = os.path.join(current_path, name)
                if not os.path.exists(folder_path):
                    os.makedirs(folder_path, exist_ok=True)
                    print(f"Folder created: {folder_path}")
                process_folder(folder_path, content)  # Recur into the subfolder
            else:  # 파일 처리
                file_path = os.path.join(current_path, name)
                with open(file_path, 'w') as f:
                    f.write(content)
                print(f"File created: {file_path}")

    # Start processing the reference folder structure
    process_folder(script_folder, reference_folders)

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
