import json
from typing import Any, Dict, List
import os
import subprocess
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
    output_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
    input_path = "user/{user_name}/data"

    # Find all outputs across all rules to determine unique inputs not present in outputs
    all_outputs = {out for rule in rules_data.values() for out in rule['output']}
    all_inputs = {inp for rule in rules_data.values() for inp in rule['input']}
    unique_inputs = [inp for inp in all_inputs if inp not in all_outputs]  # Non-duplicated inputs compared to outputs

    # Iterate through each rule in the dictionary
    for rule_id, rule in rules_data.items():
        # Start rule block with rule name
        snakemake_code += f"rule {rule['name']}:\n"
        
        # Input section
        if 'input' in rule and rule['input']:
            input_files = ",\n        ".join([
                f"{os.path.splitext(os.path.basename(inp))[0]}=\"{input_path}/{{{inp}}}\""
                for inp in rule['input']
            ])
            snakemake_code += f"    input:\n        {input_files}\n"
        
        # Output section
        if 'output' in rule and rule['output']:
            output_files = ",\n        ".join([
                f"{os.path.splitext(os.path.basename(out))[0]}=\"{output_path}/{out}\""
                for out in rule['output']
            ])
            snakemake_code += f"    output:\n        {output_files}\n"
        
        # Shell section based on script type
        if 'script' in rule and rule['script']:
            script_path = "/app/plugin/{plugin_name}/scripts/" + rule['script']
            if rule['script'].endswith('.py'):
                shell_command = f"/opt/conda/envs/snakemake/bin/python {script_path}"
            elif rule['script'].endswith('.R'):
                shell_command = f"/Rscript {script_path}"
            else:
                shell_command = script_path  # For other script types (e.g., .sh)

            # Add parameters to shell command
            if 'parameters' in rule and rule['parameters']:
                param_list = " ".join([f"{{{param['name']}}}" for param in rule['parameters']])
                shell_command = f"{shell_command} {param_list}"

            snakemake_code += f"    shell:\n        \"{shell_command}\"\n"

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

    if dependency_file_name == "requirements.txt":
        print(f"Installing dependencies from {dependency_file_name} using pip...")
        # Install dependencies using pip in the conda environment
        subprocess.run([pip_executable, "install", "-r", dependency_file_name], check=True)

    elif dependency_file_name == "environment.yml":
        print(f"Installing dependencies from {dependency_file_name} using conda...")
        # Install dependencies using conda in the conda environment
        subprocess.run([conda_executable, "env", "update", "--file", dependency_file_name, "--prune"], check=True)

    elif dependency_file_name == "renv.lock":
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


def create_dependency_folder(dependency_folder: str):
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