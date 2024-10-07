import json
from typing import Any, Dict, List
import os
import subprocess
from fastapi import HTTPException
import yaml
import re

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
    pos_x_visualize = pos_x_resultfile + pos_x_increment

    pos_y_inputfile = pos_y_start
    pos_y_datatable = pos_y_start
    pos_y_scatterplot = pos_y_start
    pos_y_algorithm = pos_y_start
    pos_y_resultfile = pos_y_start
    pos_y_visualize = pos_y_start

    # 새로운 노드를 저장할 리스트 초기화
    new_nodes = []
    resultfile_nodes = []
    visualize_nodes = []

    # 기존 노드를 순회하며 새로운 노드를 생성합니다.
    for key, node in original_data.items():
        # inputs를 확인하여 connections가 비어있는 경우 새로운 노드 생성
        for input_key, input_val in node['inputs'].items():
            if not input_val['connections']:
                index = int(input_key.split('_')[1]) - 1
                input_file = node['data']['inputs'][index]
                if input_file.endswith('.csv') or input_file.endswith('.h5ad'):
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
                    visualize_node = {
                        "id": node_id,
                        "name": "Visualize",
                        "data": {"title": output_file},
                        "class": "Visualize",
                        "html": "Visualize",
                        "typenode": "vue",
                        "inputs": {
                            "input_1": {
                                "connections": []
                            }
                        },
                        "outputs": {},
                        "pos_x": pos_x_visualize,
                        "pos_y": pos_y_visualize
                    }
                    new_drawflow["drawflow"]["Home"]["data"][str(node_id)] = visualize_node
                    visualize_nodes.append(visualize_node)
                    node_id += 1
                    pos_y_visualize += pos_y_increment
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

    # Visualize 노드의 입력을 모든 ResultFile 노드의 출력에 연결
    for visualize_node in visualize_nodes:
        visualize_node["inputs"]["input_1"]["connections"] = [
            {"node": str(resultfile_node['id']), "input": "output_1"}
            for resultfile_node in resultfile_nodes
        ]
        for resultfile_node in resultfile_nodes:
            resultfile_node["outputs"]["output_1"]["connections"].append({
                "node": str(visualize_node['id']),
                "input": "input_1"
            })

    # 생성된 new_drawflow 데이터를 출력합니다.
    print(json.dumps(new_drawflow, indent=2))
    return new_drawflow

def generate_snakemake_code(rules_data, output_file_path):
    snakemake_code = ""
    input_output_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
    logs_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/logs"
    unique_input_path = "user/{user_name}/data"

    # Find all outputs across all rules to determine unique inputs not present in outputs
    all_outputs = {out for rule in rules_data.values() for out in rule['output']}
    all_inputs = {inp for rule in rules_data.values() for inp in rule['input']}
    unique_inputs = [inp for inp in all_inputs if inp not in all_outputs]  # Non-duplicated inputs compared to outputs
    scripts_log_word = "> {log.stdout} 2> {log.stderr}"

    # Iterate through each rule in the dictionary
    for rule_id, rule in rules_data.items():
        # Start rule block with rule name
        snakemake_code += f"rule {rule['name']}:\n"
        
        # Input section
        if 'input' in rule and rule['input']:
            input_files = ",\n        ".join([
                f"{os.path.splitext(os.path.basename(inp))[0]}=\"{unique_input_path}/{{{inp}}}\"" if inp in unique_inputs
                else f"{os.path.splitext(os.path.basename(inp))[0]}=\"{input_output_path}/{inp}\""
                for inp in rule['input']
            ])
            snakemake_code += f"    input:\n        {input_files}\n"
        
        # Output section
        if 'output' in rule and rule['output']:
            output_files = ",\n        ".join([
                f"{os.path.splitext(os.path.basename(out))[0]}=\"{input_output_path}/{out}\""
                for out in rule['output']
            ])
            snakemake_code += f"    output:\n        {output_files}\n"

        # Params section
        if 'parameters' in rule and rule['parameters']:
            param_list = []
            for param in rule['parameters']:
                if param['name'] == "clusters" and param['type'] == "h5adParameter":
                    # Add clusters using ";".join for a list
                    param_list.append(f'clusters=lambda wc: ";".join({{{param["name"]}}})')
                elif param['type'] != 'inputFile' and param['type'] != 'outputFile':
                    param_list.append(f"{param['name']}={{{param['name']}}}")
            
            if param_list:
                param_list_str = ",\n        ".join(param_list)
                snakemake_code += f"    params:\n        {param_list_str}\n"

        # Log section
        snakemake_code += f"    log:\n"
        snakemake_code += f"        stdout=\"{logs_path}/{rule['name']}.stdout\",\n"
        snakemake_code += f"        stderr=\"{logs_path}/{rule['name']}.stderr\"\n"
        
        # Shell section based on script type
        if 'script' in rule and rule['script']:
            script_path = "/app/plugin/{plugin_name}/scripts/" + rule['script']
            if rule['script'].endswith('.py'):
                shell_command = f"/opt/conda/envs/snakemake/bin/python {script_path}"
            elif rule['script'].endswith('.R'):
                shell_command = f"/Rscript {script_path}"
            else:
                shell_command = script_path  # For other script types (e.g., .sh)

            # Add parameters to shell command with input/output/params distinction
            if 'parameters' in rule and rule['parameters']:
                param_list = []
                for param in rule['parameters']:
                    if param['type'] == 'inputFile':
                        param_list.append(f"{{input.{param['name']}}}")
                    elif param['type'] == 'outputFile':
                        param_list.append(f"{{output.{param['name']}}}")
                    elif param['name'] == "clusters" and param['type'] == "h5adParameter":
                        param_list.append(f"\'{{params.{param['name']}}}\'")
                    else:
                        param_list.append(f"{{params.{param['name']}}}")
                
                # Join all the parameters with a space
                param_list_str = " ".join(param_list)
                shell_command = f"{shell_command} {param_list_str}"

            snakemake_code += f"    shell:\n        \"{shell_command} {scripts_log_word}\"\n"

        # Add a newline for separation between rules
        snakemake_code += "\n"

    # Write the generated Snakemake code to a file
    with open(output_file_path, 'w') as file:
        file.write(snakemake_code)

    print(f"Snakemake code has been generated and saved to {output_file_path}.")

def install_dependencies(dependency_file_name: str):
    """
    Install dependencies based on the given file name in the /opt/conda/envs/snakemake environment.
    
    Parameters:
        dependency_file_name (str): The name of the dependency file (requirements.txt, environment.yml, renv.lock)
    """
    # Set the path to the conda environment
    conda_env_path = "/opt/conda/envs/snakemake"
    python_executable = f"{conda_env_path}/bin/python"
    pip_executable = f"{conda_env_path}/bin/pip"
    conda_executable = f"{conda_env_path}/bin/conda"
    rscript_executable = f"{conda_env_path}/bin/Rscript"

    if "requirements.txt" in dependency_file_name:
        print(f"Installing dependencies from {dependency_file_name} using pip...")
        # Install dependencies using pip in the conda environment
        subprocess.run([pip_executable, "install", "-r", dependency_file_name], check=True)

    elif "environment.yml" in dependency_file_name or "environment.yaml" in dependency_file_name:
        print(f"Installing dependencies from {dependency_file_name} using conda...")
        # Install dependencies using conda in the conda environment
        subprocess.run([conda_executable, "env", "update", "--file", dependency_file_name, "--prune"], check=True)

    elif "renv.lock" in dependency_file_name:
        print(f"Installing R dependencies from {dependency_file_name} using renv...")
        # Install R dependencies using renv in the conda environment
        subprocess.run([rscript_executable, "-e", f'renv::restore(lockfile = "{dependency_file_name}")'], check=True)

    else:
        print(f"Unsupported file: {dependency_file_name}. Please provide a valid dependency file.")
        return

    print(f"Dependencies from {dependency_file_name} have been installed successfully.")

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

def normalize_pkg_name(name: str) -> str:
    """Normalize package names for consistent comparison."""
    return name.lower().replace("-", "_")

def check_requirements_txt(requirements_file: str):
    """Check if dependencies in requirements.txt are installed."""
    
    # Check if the requirements file exists
    if not os.path.exists(requirements_file):
        raise FileNotFoundError(f"{requirements_file} not found.")
    
    # Get installed packages via pip freeze
    installed_packages = subprocess.run(
        ["/opt/conda/envs/snakemake/bin/pip", "freeze"],
        stdout=subprocess.PIPE,
        text=True
    ).stdout.splitlines()

    # Read the requirements.txt file
    with open(requirements_file, 'r') as f:
        required_packages = f.read().splitlines()

    # Parse installed packages to a dict
    installed_dict = {normalize_pkg_name(re.split(r"==|>=|<=|~=", pkg)[0]): pkg for pkg in installed_packages}

    installed = []
    missing = []

    # Check each requirement
    for req in required_packages:
        pkg_name = normalize_pkg_name(re.split(r"==|>=|<=|~=", req)[0])
        if pkg_name in installed_dict:
            installed.append(pkg_name)
        else:
            missing.append(pkg_name)

    # Return result in dict form
    if len(missing) > 0:
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
    
    # Get installed conda packages via conda list
    installed_conda_packages = subprocess.run(
        ["/opt/conda/envs/snakemake/bin/conda", "list"],
        stdout=subprocess.PIPE,
        text=True
    ).stdout.splitlines()

    # Get installed pip packages via pip freeze
    installed_pip_packages = subprocess.run(
        ["/opt/conda/envs/snakemake/bin/pip", "freeze"],
        stdout=subprocess.PIPE,
        text=True
    ).stdout.splitlines()

    # Parse installed pip packages to a dict
    installed_pip_dict = {pkg.split("==")[0] for pkg in installed_pip_packages}

    # Parse the environment.yml file using the yaml module
    try:
        with open(environment_file, 'r') as f:
            env_data = yaml.safe_load(f)  # Parse the environment.yml
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing {environment_file}: {e}")

    required_conda_packages = []
    required_pip_packages = []

    # Ensure there are dependencies to process
    if 'dependencies' in env_data:
        for dep in env_data['dependencies']:
            # Check for pip section, if it's a pip section add to pip package list
            if isinstance(dep, dict) and 'pip' in dep:
                # Ensure that 'dep['pip']' is a list of strings, handle each package individually
                for pip_pkg in dep['pip']:
                    if isinstance(pip_pkg, str):
                        required_pip_packages.append(pip_pkg)
            # Otherwise, it's a conda package
            elif isinstance(dep, str):
                required_conda_packages.append(normalize_pkg_name(dep.split("=")[0].strip()))

    installed_dict = {normalize_pkg_name(pkg.split()[0]) for pkg in installed_conda_packages if len(pkg.split()) > 0}

    conda_installed = []
    conda_missing = []
    pip_installed = []
    pip_missing = []

    # Check Conda packages
    for req in required_conda_packages:
        pkg_name = normalize_pkg_name(req.split("=")[0].strip())
        if pkg_name in installed_dict:
            conda_installed.append(pkg_name)
        else:
            conda_missing.append(pkg_name)

    # Check Pip packages
    for req in required_pip_packages:
        pkg_name = normalize_pkg_name(req.split("==")[0].strip())
        if pkg_name in installed_pip_dict:
            pip_installed.append(pkg_name)
        else:
            pip_missing.append(pkg_name)

    # Combine missing packages
    missing_packages = conda_missing + pip_missing

    # Return result in dict form
    if len(missing_packages) > 0:
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
    else:
        print(f"Unsupported dependency file: {dependency_file_name}")
        return {
            "installed_status": False,
            "message": f"Unsupported file: {dependency_file_name}"
        }
