import os
import pandas as pd
import re
import scanpy as sc
import numpy as np
import logging
import json
import gc
from typing import Dict, List, Optional, Tuple, Any
from app.core.exceptions import SnakefileGenerationError

# Configure logging
logger = logging.getLogger(__name__)

# Custom exceptions for better error handling
class WorkflowUtilsError(Exception):
    """Base exception for workflow utilities."""
    pass

class FileValidationError(WorkflowUtilsError):
    """Exception raised when file validation fails."""
    def __init__(self, message: str, file_path: str = None, context: Dict[str, Any] = None):
        super().__init__(message)
        self.file_path = file_path
        self.context = context or {}
        logger.error(f"File validation error: {message}. Path: {file_path}. Context: {context}")

class SnakefileProcessingError(WorkflowUtilsError):
    """Exception raised when Snakefile processing fails."""
    def __init__(self, message: str, snakefile_path: str = None, rule_name: str = None):
        super().__init__(message)
        self.snakefile_path = snakefile_path
        self.rule_name = rule_name
        logger.error(f"Snakefile processing error: {message}. File: {snakefile_path}, Rule: {rule_name}")

class AlgorithmResolutionError(WorkflowUtilsError):
    """Exception raised when algorithm path resolution fails."""
    def __init__(self, message: str, node_id: str = None, available_files: List = None):
        super().__init__(message)
        self.node_id = node_id
        self.available_files = available_files
        logger.error(f"Algorithm resolution error: {message}. Node ID: {node_id}")

# def load_tab_file(file_path: str, max_rows: int = 5000):
def load_tab_file(file_path: str):
    """
    Load tabular data from CSV, TSV, or H5AD files.
    
    Args:
        file_path: Path to the file to load
        
    Returns:
        pandas.DataFrame: Loaded data
        
    Raises:
        FileValidationError: If file doesn't exist or has invalid extension
        WorkflowUtilsError: If file loading fails
    """
    logger.info(f"Loading tabular file: {file_path}")
    
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileValidationError(
            f"File does not exist",
            file_path=file_path,
            context={"operation": "load_tab_file"}
        )
    
    # Extract the file extension
    file_extension = os.path.splitext(file_path)[1].lower()
    
    # Check if the file extension is supported
    supported_extensions = ['.csv', '.tsv', '.h5ad']
    if file_extension not in supported_extensions:
        raise FileValidationError(
            f"Unsupported file extension. Supported: {', '.join(supported_extensions)}",
            file_path=file_path,
            context={"extension": file_extension, "supported": supported_extensions}
        )
    
    try:
        if file_extension == '.h5ad':
            logger.debug(f"Loading H5AD file: {file_path}")
            # 메모리 누수 방지를 위해 adata를 명시적으로 해제
            adata = None
            try:
                adata = sc.read_h5ad(file_path)
                # 기본 데이터프레임으로 adata.obs 사용 (.copy()로 참조 끊기)
                df = adata.obs.copy()

                # X_umap이 있는 경우에만 추가
                if 'X_umap' in adata.obsm:
                    umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)
                    df = pd.concat([df, umap_df], axis=1)
                # X_pca가 있는 경우 대체 사용
                elif 'X_pca' in adata.obsm:
                    pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :2], columns=['X', 'Y'], index=adata.obs.index)
                    df = pd.concat([df, pca_df], axis=1)
                # 둘 다 없는 경우 raw 데이터에서 처음 두 컬럼 사용
                else:
                    # sparse matrix인 경우 dense로 변환
                    from scipy.sparse import issparse
                    # 데이터가 1차원인 경우 처리
                    if adata.X.shape[1] == 1:
                        col_data = adata.X[:, 0]
                        if issparse(col_data):
                            col_data = col_data.toarray().ravel()
                        raw_df = pd.DataFrame({
                            'X': col_data,
                            'Y': np.zeros(adata.X.shape[0])  # Y값을 0으로 설정
                        }, index=adata.obs.index)
                    else:
                        col_data = adata.X[:, :2]
                        if issparse(col_data):
                            col_data = col_data.toarray()
                        raw_df = pd.DataFrame(col_data, columns=['X', 'Y'], index=adata.obs.index)
                    df = pd.concat([df, raw_df], axis=1)

                logger.debug(f"Successfully loaded H5AD file with {df.shape[0]} observations")
                return df
            finally:
                if adata is not None:
                    del adata
                gc.collect()
        else:
            # Determine the separator based on the file extension
            sep = ',' if file_extension == '.csv' else '\t'
            logger.debug(f"Loading {'CSV' if file_extension == '.csv' else 'TSV'} file: {file_path}")
            
            # Load the file with pandas
            df = pd.read_csv(file_path, sep=sep)
            logger.debug(f"Successfully loaded tabular file with shape {df.shape}")
            
            # Return the dataframe
            return df
    except Exception as e:
        logger.error(f"Failed to load file {file_path}: {str(e)}")
        raise WorkflowUtilsError(f"Failed to load file: {str(e)}") from e

def transform_df_to_vgt_format(df):
    columns = [
        {
            "label": col,
            "field": col,
            "type": "number" if df[col].dtype in [int, float] else "string"
        }
        for col in df.columns
    ]

    rows = df.to_dict(orient='records')
    
    return {"columns": columns, "rows": rows}

def extract_all_algorithms(workflow_info):
    algorithms = []
    
    # 탐색하여 class가 "Algorithm"인 모든 객체를 찾음
    for key, value in workflow_info.items():
        if value.get("class") == "Algorithm":
            algorithm_data = value.get("data", {})
            algorithm_id = value.get("id")

            algorithms.append({
                "id": algorithm_id,
                "files": algorithm_data.get("files"),
                "selectedPlugin": algorithm_data.get("selectedPlugin"),
                "selectedPluginInputOutput": algorithm_data.get("selectedPluginInputOutput"),
                "selectedPluginRules": algorithm_data.get("selectedPluginRules")
            })

    return algorithms


def find_connected_visualization_nodes(workflow_data: dict, algorithm_id) -> List[str]:
    """
    Find all Visualization node IDs connected to a given Algorithm node.

    Traverses workflow graph from Algorithm node through ResultFiles to find
    Visualization nodes.

    Args:
        workflow_data: The workflow data dictionary (drawflow.Home.data structure)
        algorithm_id: The ID of the Algorithm node to start from

    Returns:
        List of Visualization node IDs connected to the Algorithm node
    """
    visualization_node_ids = []
    algorithm_id_str = str(algorithm_id)

    # Find the Algorithm node
    algorithm_node = None
    for key, node in workflow_data.items():
        if str(node.get("id")) == algorithm_id_str and node.get("class") == "Algorithm":
            algorithm_node = node
            break

    if not algorithm_node:
        logger.warning(f"Algorithm node with ID {algorithm_id} not found")
        return visualization_node_ids

    # Get direct connections from Algorithm node outputs
    outputs = algorithm_node.get("outputs", {})
    connected_node_ids = []

    for output_key, output_data in outputs.items():
        connections = output_data.get("connections", [])
        for connection in connections:
            connected_node_id = connection.get("node")
            if connected_node_id:
                connected_node_ids.append(str(connected_node_id))

    # Process connected nodes
    for connected_id in connected_node_ids:
        connected_node = workflow_data.get(connected_id)
        if not connected_node:
            continue

        node_class = connected_node.get("class", "")

        if node_class == "Visualization":
            visualization_node_ids.append(connected_id)
        elif node_class == "ResultFiles":
            # Follow ResultFiles outputs to find Visualization nodes
            rf_outputs = connected_node.get("outputs", {})
            for rf_output_key, rf_output_data in rf_outputs.items():
                rf_connections = rf_output_data.get("connections", [])
                for rf_connection in rf_connections:
                    rf_connected_id = rf_connection.get("node")
                    if rf_connected_id:
                        rf_connected_node = workflow_data.get(str(rf_connected_id))
                        if rf_connected_node and rf_connected_node.get("class") == "Visualization":
                            visualization_node_ids.append(str(rf_connected_id))

    # Remove duplicates while preserving order
    seen = set()
    unique_ids = []
    for vid in visualization_node_ids:
        if vid not in seen:
            seen.add(vid)
            unique_ids.append(vid)

    logger.debug(f"Found {len(unique_ids)} Visualization nodes connected to Algorithm {algorithm_id}")
    return unique_ids


def extract_algorithm_data(workflow_info, algorithm_id):
    # 탐색하여 class가 "Algorithm"이고 id가 algorithm_id인 객체를 찾음
    for key, value in workflow_info.items():
        if value.get("class") == "Algorithm" and str(value.get("id")) == algorithm_id:
            algorithm_data = value.get("data", {})
            
            # 필드 추출
            selected_plugin = algorithm_data.get("selectedPlugin")
            selected_plugin_input_output = algorithm_data.get("selectedPluginInputOutput")
            selected_plugin_rules = algorithm_data.get("selectedPluginRules")

            # 결과 반환
            return {
                "id": algorithm_id,
                "files": algorithm_data.get("files"),
                "selectedPlugin": selected_plugin,
                "selectedPluginInputOutput": selected_plugin_input_output,
                "selectedPluginRules": selected_plugin_rules
            }
    
    # "Algorithm" 클래스를 가진 객체가 없을 경우 빈 딕셔너리 반환
    return {}

def extract_visualization_data(workflow_info, node_id):
    # 탐색하여 class가 "Visualization"이고 id가 node_id인 객체를 찾음
    for key, value in workflow_info.items():
        print(value.get("class"))
        print(value.get("id"))
        if value.get("class") == "Visualization" and str(value.get("id")) == node_id:
            visualization_data = value.get("data", {})
            
            # 필드 추출
            selected_visualization_params= visualization_data.get("selectedVisualizationParams")
            selected_visualization_title = visualization_data.get("selectedVisualizationTitle")

            # 결과 반환
            return {
                "selectedVisualizationParams": selected_visualization_params,
                "selectedVisualizationTitle": selected_visualization_title
            }
    
    # "Visualization" 클래스를 가진 객체가 없을 경우 빈 딕셔너리 반환
    return {}

def generate_user_input(selectedPluginInputOutput, username=None, file_sources=None):
    """
    플러그인 입출력 파라미터에서 사용자 입력 파일 매핑을 생성한다.

    Args:
        selectedPluginInputOutput: 플러그인 입출력 파라미터 배열
        username: 유저명 (전체 경로 구성 시 사용)
        file_sources: {파라미터키: "user"|"shared"} 매핑 (None이면 기존 동작)

    Returns:
        {"input.h5ad": "user/john/data/pbmc.h5ad"} 또는
        {"input.h5ad": "tutorials/pbmc.h5ad"} 또는
        {"input.h5ad": "pbmc.h5ad"} (하위 호환: username/file_sources 미제공 시)
    """
    user_input = {}

    for parameter in selectedPluginInputOutput:
        if parameter.get("type") in ("inputFile", "optionalInputFile"):
            parameter_key = parameter.get("defaultValue")
            file_name = parameter.get("file_name")

            if not file_name:
                continue

            # 하위 호환: username이 없으면 기존 동작 (파일명만 반환)
            if username is None or file_sources is None:
                user_input[parameter_key] = file_name
            else:
                source = file_sources.get(parameter_key, "user")
                if source == "shared":
                    user_input[parameter_key] = f"tutorials/{file_name}"
                else:
                    user_input[parameter_key] = f"user/{username}/data/{file_name}"

    return user_input


def extract_file_sources(algorithm_data):
    """Algorithm 노드에서 file_source 매핑 추출.

    selectedPluginInputOutput 배열의 각 inputFile 파라미터에서
    fileSource 필드를 읽어 {defaultValue: source} 매핑을 반환한다.

    Args:
        algorithm_data: algorithm dict (selectedPluginInputOutput 포함)

    Returns:
        {"input.h5ad": "shared", "geneList.txt": "user", ...}
    """
    sources = {}
    input_output = algorithm_data.get('selectedPluginInputOutput', [])
    for param in input_output:
        if param.get('type') in ('inputFile', 'optionalInputFile'):
            key = param.get('defaultValue')
            source = param.get('fileSource', 'user')
            sources[key] = source
    return sources

def generate_plugin_params(selectedPluginRules):
    plugin_params = {}

    # selectedPluginRules에서 파라미터 추출 및 사용자 입력에 추가
    for rule in selectedPluginRules:
        for parameter in rule.get("parameters", []):
            # 파라미터 이름 추출
            parameter_name = parameter.get("name")
            
            plugin_params[parameter_name] = parameter.get("defaultValue")

    return plugin_params

def generate_visualization_params(selectedVisualizationParams):
    visualization_params = {}
    visualization_inputs = {}
    visualization_outputs = {}

    # 모든 파라미터를 적절한 딕셔너리에 분류
    for parameter in selectedVisualizationParams:
        param_type = parameter.get("type")
        param_name = parameter.get("name")
        param_value = parameter.get("defaultValue")
        
        # inputFile과 optionalInputFile은 visualization_inputs에 할당
        if param_type in ["inputFile", "optionalInputFile"]:
            if "target" in param_name:
                param_name = param_name + parameter.get("fileExtension")
            visualization_inputs[param_name] = param_value
            
        # outputFile은 visualization_outputs에 할당
        elif param_type == "outputFile":
            visualization_outputs[param_name] = param_value
            
        # 나머지는 visualization_params에 할당
        else:
            visualization_params[param_name] = param_value

    return visualization_params, visualization_inputs, visualization_outputs

def extract_target_data(selectedPluginInputOutput, user_workflow_task_path):
    data_list = []

    user_workflow_task_path = os.path.join(user_workflow_task_path, "results")

    # type이 "output"인 데이터들 중, activated가 True인 데이터들을 찾아서 리스트에 추가
    for data in selectedPluginInputOutput:
        if data.get("type") == "outputFile" and data.get("activate"):
            data_list.append(os.path.join(user_workflow_task_path, data.get("defaultValue")))

    # 데이터 리스트 반환
    return data_list

def extract_rule_block(snakefile_path: str, visualization_title: str, result_path: str, 
                      parameters: Optional[Dict[str, Any]] = None, 
                      file_mappings: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """
    Enhanced rule extraction with dynamic parameter replacement and better error handling.
    
    Args:
        snakefile_path: Path to the Snakefile containing rules
        visualization_title: Name of the rule to extract
        result_path: Path where extracted rule should be saved
        parameters: Optional dict of parameters to replace in the rule
        file_mappings: Optional dict mapping file placeholders to actual paths
        
    Returns:
        Tuple of (extracted_content, result_path)
        
    Raises:
        SnakefileProcessingError: If snakefile processing fails
        FileValidationError: If files don't exist or are invalid
    """
    logger.info(f"Extracting rule '{visualization_title}' from {snakefile_path}")
    
    # Validate inputs
    if not snakefile_path:
        raise SnakefileProcessingError(
            "Snakefile path is required",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )
    
    if not visualization_title:
        raise SnakefileProcessingError(
            "Rule name is required",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )
    
    if not os.path.exists(snakefile_path):
        raise FileValidationError(
            f"Snakefile not found",
            file_path=snakefile_path,
            context={"operation": "extract_rule_block", "rule": visualization_title}
        )
    
    if not os.path.isfile(snakefile_path):
        raise FileValidationError(
            f"Path is not a file",
            file_path=snakefile_path,
            context={"operation": "extract_rule_block", "rule": visualization_title}
        )
    
    # rule 블록들을 분리하기 위한 시작 패턴
    rule_pattern = f"rule {visualization_title}:"
    
    # 파일 내용을 줄 단위로 분리
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except UnicodeDecodeError as e:
        raise FileValidationError(
            f"Failed to read Snakefile: encoding error",
            file_path=snakefile_path,
            context={"error": str(e), "encoding": "utf-8"}
        )
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to read Snakefile: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        ) from e
    
    result = []
    is_target_rule = False
    current_indent = 0
    rule_found = False
    
    for line_num, line in enumerate(lines, 1):
        # 대상 rule 블록의 시작을 찾음
        if line.startswith(rule_pattern):
            is_target_rule = True
            rule_found = True
            current_indent = len(line) - len(line.lstrip())
            result.append(line)
            logger.debug(f"Found rule '{visualization_title}' at line {line_num}")
            continue
            
        # rule 블록 내부의 내용을 수집
        if is_target_rule:
            # 현재 줄의 들여쓰기 수준 확인
            if line.strip() == "":
                result.append(line)  # Preserve empty lines within the rule
                continue
                
            line_indent = len(line) - len(line.lstrip())
            
            # 들여쓰기가 같거나 큰 경우에만 해당 rule에 속한 것으로 간주
            if line_indent > current_indent:
                result.append(line)
            else:
                break
    
    if not rule_found:
        try:
            available_rules = _get_available_rules(snakefile_path)
        except Exception:
            available_rules = ["Unable to determine available rules"]
        raise SnakefileProcessingError(
            f"Rule '{visualization_title}' not found. Available rules: {available_rules}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )
                
    # 빈 줄 제거하고 결과 문자열 생성
    result = [line.rstrip() for line in result if line.strip() or line.strip() == ""]
    extracted_content = 'import os\n' + '\n'.join(result)
    
    # Apply parameter and file mapping replacements if provided
    if parameters or file_mappings:
        extracted_content = _apply_replacements(extracted_content, parameters, file_mappings)
    
    # 추출한 내용을 파일에 쓰기
    try:
        # Ensure parent directory exists
        result_dir = os.path.dirname(result_path)
        if result_dir:
            os.makedirs(result_dir, exist_ok=True)
        
        with open(result_path, 'w', encoding='utf-8') as file:
            file.write(extracted_content)
        logger.info(f"Successfully extracted rule to {result_path}")
    except PermissionError as e:
        raise FileValidationError(
            f"Permission denied writing to result path",
            file_path=result_path,
            context={"operation": "write_extracted_rule", "error": str(e)}
        )
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to write extracted rule: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        ) from e

    return extracted_content, result_path

def resolve_algorithm_path_from_files(workflow_info: dict, current_node_id: str, available_files: list) -> Optional[str]:
    """
    Resolve Algorithm node path from visualization node connections and available files.
    
    This function traverses the workflow graph to find the Algorithm node connected to a 
    Visualization node through ResultFiles nodes. The connection flow is typically:
    Algorithm -> ResultFiles -> Visualization
    
    Args:
        workflow_info (dict): Complete workflow information containing nodes and connections
                             Expected structure: {node_key: {id, class, inputs, outputs, ...}}
        current_node_id (str): ID of the visualization node requesting files
        available_files (list): List of available files from ResultFiles node
                               Can be list of strings or list of dicts with 'name'/'filename' keys
        
    Returns:
        Optional[str]: Algorithm ID that can be used to construct file paths like 
                      'algorithm_{id}/results/{filename}', or None if not found
        
    Raises:
        AlgorithmResolutionError: If algorithm path resolution fails
        
    Example:
        algorithm_id = resolve_algorithm_path_from_files(
            workflow_info=workflow.workflow_info,
            current_node_id="3",
            available_files=[{"name": "result.csv"}, {"name": "output.txt"}]
        )
        # Returns: "1" (if algorithm node has id=1)
    """
    try:
        # Validate inputs
        if not workflow_info:
            raise AlgorithmResolutionError(
                "Workflow info is required",
                node_id=current_node_id
            )
        
        if not current_node_id:
            raise AlgorithmResolutionError(
                "Current node ID is required",
                available_files=available_files
            )
        
        logger.info(f"Resolving algorithm path for visualization node {current_node_id}")
        
        # Find the visualization node
        visualization_node = None
        for node_key, node_value in workflow_info.items():
            if str(node_value.get("id")) == str(current_node_id):
                if node_value.get("class") == "Visualization" or node_value.get("data", {}).get("isVisualization", False):
                    visualization_node = node_value
                    break
        
        if not visualization_node:
            raise AlgorithmResolutionError(
                f"Visualization node with ID {current_node_id} not found in workflow",
                node_id=current_node_id,
                available_files=available_files
            )
        
        logger.debug(f"Found visualization node: {visualization_node.get('html', 'Unknown')}")
        
        # Get input connections from the visualization node
        inputs = visualization_node.get("inputs", {})
        if not inputs:
            logger.warning(f"No inputs found for visualization node {current_node_id}")
            return None
            
        # Traverse input connections to find ResultFiles nodes
        resultfiles_node_ids = []
        for input_key, input_data in inputs.items():
            connections = input_data.get("connections", [])
            for connection in connections:
                connected_node_id = connection.get("node")
                if connected_node_id:
                    # Find the connected node in workflow_info
                    for node_key, node_value in workflow_info.items():
                        if str(node_value.get("id")) == str(connected_node_id):
                            # Check if this is a ResultFiles node
                            node_html = node_value.get("html", "")
                            if "ResultFiles" in node_html or "result" in node_html.lower():
                                resultfiles_node_ids.append(connected_node_id)
                                logger.debug(f"Found ResultFiles node: {connected_node_id}")
                            break
        
        if not resultfiles_node_ids:
            logger.warning(f"No ResultFiles nodes connected to visualization node {current_node_id}")
            return None
            
        # From ResultFiles nodes, find connected Algorithm nodes
        algorithm_ids = []
        for resultfiles_id in resultfiles_node_ids:
            for node_key, node_value in workflow_info.items():
                if str(node_value.get("id")) == str(resultfiles_id):
                    # Get input connections from ResultFiles node
                    rf_inputs = node_value.get("inputs", {})
                    for input_key, input_data in rf_inputs.items():
                        connections = input_data.get("connections", [])
                        for connection in connections:
                            algorithm_node_id = connection.get("node")
                            if algorithm_node_id:
                                # Verify this is an Algorithm node
                                for alg_key, alg_value in workflow_info.items():
                                    if str(alg_value.get("id")) == str(algorithm_node_id):
                                        if alg_value.get("class") == "Algorithm":
                                            algorithm_ids.append(algorithm_node_id)
                                            logger.debug(f"Found Algorithm node: {algorithm_node_id}")
                                        break
                    break
        
        if not algorithm_ids:
            logger.warning("No Algorithm nodes found connected to ResultFiles nodes")
            return None
            
        # Try each algorithm ID to see which one has valid file paths
        for algorithm_id in algorithm_ids:
            logger.debug(f"Checking algorithm {algorithm_id} for valid file paths")
            if available_files:
                # Extract just filenames from available_files (they may include full paths)
                filenames = []
                for file_item in available_files:
                    if isinstance(file_item, dict):
                        filename = file_item.get("name") or file_item.get("filename")
                    else:
                        filename = file_item
                    
                    if filename:
                        # Extract just the filename without path
                        filenames.append(os.path.basename(filename))
                
                if filenames:
                    try:
                        # This is a placeholder check - we'll validate paths in validate_file_paths
                        logger.info(f"Successfully resolved algorithm path: {algorithm_id}")
                        return str(algorithm_id)
                    except Exception as e:
                        logger.debug(f"Algorithm {algorithm_id} validation failed: {e}")
                        continue
        
        raise AlgorithmResolutionError(
            "No valid algorithm path found with existing files",
            node_id=current_node_id,
            available_files=available_files
        )
        
    except AlgorithmResolutionError:
        # Re-raise our custom exceptions
        raise
    except Exception as e:
        raise AlgorithmResolutionError(
            f"Unexpected error resolving algorithm path: {str(e)}",
            node_id=current_node_id,
            available_files=available_files
        ) from e

def validate_file_paths(user_path: str, workflow_id: int, algorithm_id: str, filenames: list) -> list:
    """
    Validate that files exist in the expected algorithm results directory.
    
    This function constructs the expected file paths following the CellCraft directory structure:
    {user_path}/workflow_{workflow_id}/algorithm_{algorithm_id}/results/{filename}
    
    Args:
        user_path (str): Base user data path (e.g., backend/user_data/{user_id}/)
        workflow_id (int): Workflow ID for path construction
        algorithm_id (str): Algorithm ID for path construction  
        filenames (list): List of filenames to validate (just filenames, not full paths)
        
    Returns:
        list: List of validated absolute file paths that exist on the filesystem
        
    Raises:
        FileValidationError: If validation fails or files don't exist
        
    Example:
        validated_paths = validate_file_paths(
            user_path="/home/user/cellcraft/backend/user_data/123/",
            workflow_id=5,
            algorithm_id="2", 
            filenames=["result.csv", "output.txt"]
        )
        # Returns: ["/home/user/.../workflow_5/algorithm_2/results/result.csv", ...]
    """
    try:
        # Validate inputs
        if not user_path:
            raise FileValidationError(
                "User path is required",
                context={"workflow_id": workflow_id, "algorithm_id": algorithm_id}
            )
        
        if not workflow_id:
            raise FileValidationError(
                "Workflow ID is required",
                context={"user_path": user_path, "algorithm_id": algorithm_id}
            )
        
        if not algorithm_id:
            raise FileValidationError(
                "Algorithm ID is required",
                context={"user_path": user_path, "workflow_id": workflow_id}
            )
        
        logger.info(f"Validating file paths for algorithm {algorithm_id} in workflow {workflow_id}")
        
        if not filenames:
            logger.warning("No filenames provided for validation")
            return []
        
        # Construct the expected results directory path
        results_dir = os.path.join(user_path, f"workflow_{workflow_id}", f"algorithm_{algorithm_id}", "results")
        logger.debug(f"Expected results directory: {results_dir}")
        
        # Check if the results directory exists
        if not os.path.exists(results_dir):
            raise FileValidationError(
                f"Results directory not found",
                file_path=results_dir,
                context={
                    "workflow_id": workflow_id,
                    "algorithm_id": algorithm_id,
                    "expected_dir": results_dir
                }
            )
        
        if not os.path.isdir(results_dir):
            raise FileValidationError(
                f"Results path is not a directory",
                file_path=results_dir,
                context={
                    "workflow_id": workflow_id,
                    "algorithm_id": algorithm_id
                }
            )
        
        validated_paths = []
        missing_files = []
        
        # Validate each filename
        for filename in filenames:
            if not filename:
                continue
                
            # Extract just the filename without any path components
            clean_filename = os.path.basename(filename)
            expected_path = os.path.join(results_dir, clean_filename)
            
            if os.path.exists(expected_path) and os.path.isfile(expected_path):
                # Return relative path instead of absolute path to avoid Docker path issues
                validated_paths.append(expected_path)
                logger.debug(f"Validated file: {expected_path}")
            else:
                missing_files.append(clean_filename)
                logger.warning(f"File not found: {expected_path}")
        
        if not validated_paths:
            raise FileValidationError(
                f"No valid files found in results directory",
                file_path=results_dir,
                context={
                    "workflow_id": workflow_id,
                    "algorithm_id": algorithm_id,
                    "missing_files": missing_files,
                    "requested_files": filenames
                }
            )
        
        if missing_files:
            logger.warning(f"Some files were not found: {', '.join(missing_files)}")
        
        logger.info(f"Successfully validated {len(validated_paths)} file paths")
        return validated_paths
        
    except FileValidationError:
        # Re-raise file validation errors as-is
        raise
    except Exception as e:
        raise FileValidationError(
            f"Unexpected error during file validation: {str(e)}",
            context={
                "user_path": user_path,
                "workflow_id": workflow_id,
                "algorithm_id": algorithm_id,
                "filenames": filenames
            }
        ) from e


def generate_visualization_snakefile(plugin_path: str, selected_script_name: str, 
                                   visualization_params: Dict[str, Any], 
                                   file_mappings: Dict[str, str], 
                                   output_path: str) -> str:
    """
    Generate a standalone executable Snakefile for a specific visualization script.
    
    This function reads a plugin's Snakefile, extracts only the rule for the selected script,
    replaces parameter placeholders with actual values, and generates a standalone Snakefile
    that can be executed independently.
    
    Args:
        plugin_path: Path to the plugin directory containing Snakefile
        selected_script_name: Name of the visualization rule/script to extract
        visualization_params: Dictionary of visualization parameters and their values
        file_mappings: Dictionary mapping input file placeholders to actual file paths
        output_path: Path where the generated Snakefile should be saved
        
    Returns:
        str: Path to the generated Snakefile
        
    Raises:
        SnakefileProcessingError: If Snakefile generation fails
        FileValidationError: If plugin or files don't exist
        
    Example:
        snakefile_path = generate_visualization_snakefile(
            plugin_path="/path/to/GRNViz",
            selected_script_name="Heatmap",
            visualization_params={"Top_Genes": 20, "Sample_Size": 150},
            file_mappings={"expression.csv": "/data/expr.csv", "input.txt": "/data/network.txt"},
            output_path="/tmp/heatmap_snakefile.smk"
        )
    """
    logger.info(f"Generating visualization Snakefile for '{selected_script_name}' from {plugin_path}")
    
    # Validate inputs
    if not plugin_path:
        raise SnakefileProcessingError(
            "Plugin path is required",
            rule_name=selected_script_name
        )
    
    if not selected_script_name:
        raise SnakefileProcessingError(
            "Script name is required",
            snakefile_path=plugin_path
        )
    
    if not output_path:
        raise SnakefileProcessingError(
            "Output path is required",
            snakefile_path=plugin_path,
            rule_name=selected_script_name
        )
    
    # Validate plugin path
    if not os.path.exists(plugin_path):
        raise FileValidationError(
            f"Plugin path not found",
            file_path=plugin_path,
            context={"script_name": selected_script_name, "operation": "generate_snakefile"}
        )
    
    if not os.path.isdir(plugin_path):
        raise FileValidationError(
            f"Plugin path is not a directory",
            file_path=plugin_path,
            context={"script_name": selected_script_name}
        )
    
    # Construct Snakefile path
    snakefile_path = os.path.join(plugin_path, "Snakefile")
    if not os.path.exists(snakefile_path):
        raise FileValidationError(
            f"Snakefile not found in plugin",
            file_path=snakefile_path,
            context={"plugin_path": plugin_path, "script_name": selected_script_name}
        )
    
    # Validate the rule exists before extraction
    if not validate_visualization_rule(snakefile_path, selected_script_name):
        try:
            available_rules = _get_available_rules(snakefile_path)
        except Exception:
            available_rules = ["Unable to determine available rules"]
        raise SnakefileProcessingError(
            f"Visualization rule '{selected_script_name}' not found. Available rules: {available_rules}",
            snakefile_path=snakefile_path,
            rule_name=selected_script_name
        )
    
    # Extract the specific rule with parameter replacements
    try:
        extracted_content, _ = extract_rule_block(
            snakefile_path=snakefile_path,
            visualization_title=selected_script_name,
            result_path=output_path,
            parameters=visualization_params,
            file_mappings=file_mappings
        )
        
        logger.info(f"Successfully generated visualization Snakefile: {output_path}")
        return output_path
        
    except (SnakefileProcessingError, FileValidationError):
        # Re-raise our custom exceptions
        raise
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to generate visualization Snakefile: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=selected_script_name
        ) from e


def validate_visualization_rule(snakefile_path: str, rule_name: str) -> bool:
    """
    Check if a specific rule name exists in a Snakefile and validate its structure.
    
    Args:
        snakefile_path: Path to the Snakefile to check
        rule_name: Name of the rule to validate
        
    Returns:
        bool: True if rule exists and has valid structure, False otherwise
        
    Raises:
        FileNotFoundError: If snakefile_path doesn't exist
    """
    logger.debug(f"Validating rule '{rule_name}' in {snakefile_path}")
    
    if not os.path.exists(snakefile_path):
        raise FileNotFoundError(f"Snakefile not found: {snakefile_path}")
    
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Check if rule exists
        rule_pattern = f"rule {rule_name}:"
        if rule_pattern not in content:
            logger.debug(f"Rule '{rule_name}' not found in Snakefile")
            return False
        
        # Extract rule metadata for validation
        metadata = get_rule_metadata(snakefile_path, rule_name)
        if not metadata:
            logger.warning(f"Rule '{rule_name}' found but metadata extraction failed")
            return False
        
        # Basic structure validation
        required_sections = ['input', 'output']
        for section in required_sections:
            if section not in metadata:
                logger.warning(f"Rule '{rule_name}' missing required section: {section}")
                return False
        
        logger.debug(f"Rule '{rule_name}' validation successful")
        return True
        
    except Exception as e:
        logger.error(f"Error validating rule '{rule_name}': {e}")
        return False


def get_rule_metadata(snakefile_path: str, rule_name: str) -> Optional[Dict[str, Any]]:
    """
    Extract rule metadata including inputs, outputs, and parameters.
    
    Args:
        snakefile_path: Path to the Snakefile
        rule_name: Name of the rule to analyze
        
    Returns:
        Dict containing rule metadata or None if extraction fails
        
    Example:
        metadata = get_rule_metadata("/path/to/Snakefile", "Heatmap")
        # Returns: {
        #     "inputs": ["expression.csv", "trajectory.txt", "input.txt"],
        #     "outputs": ["Heatmap.json"],
        #     "parameters": ["Top_Genes", "Sample_Size"],
        #     "shell_command": "Rscript ..."
        # }
    """
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        rule_pattern = f"rule {rule_name}:"
        metadata = {
            'inputs': [],
            'outputs': [], 
            'parameters': [],
            'shell_command': None
        }
        
        is_target_rule = False
        current_section = None
        current_indent = 0
        
        for line in lines:
            if line.startswith(rule_pattern):
                is_target_rule = True
                current_indent = len(line) - len(line.lstrip())
                continue
            
            if not is_target_rule:
                continue
            
            if line.strip() == "":
                continue
                
            line_indent = len(line) - len(line.lstrip())
            
            # End of rule block
            if line_indent <= current_indent and line.strip():
                break
            
            stripped_line = line.strip()
            
            # Identify sections
            if stripped_line.endswith(':'):
                section_name = stripped_line[:-1]
                if section_name in ['input', 'output', 'params', 'shell']:
                    current_section = section_name
                continue
            
            # Parse section content
            if current_section == 'input':
                # Extract input file patterns
                if '=' in stripped_line:
                    file_pattern = stripped_line.split('=')[1].strip().strip('"')
                    # Extract filename from pattern
                    filename = _extract_filename_from_pattern(file_pattern)
                    if filename:
                        metadata['inputs'].append(filename)
            
            elif current_section == 'output':
                # Extract output file patterns
                if '=' in stripped_line:
                    file_pattern = stripped_line.split('=')[1].strip().strip('"')
                    filename = _extract_filename_from_pattern(file_pattern)
                    if filename:
                        metadata['outputs'].append(filename)
            
            elif current_section == 'params':
                # Extract parameter names
                if '=' in stripped_line:
                    param_name = stripped_line.split('=')[0].strip()
                    metadata['parameters'].append(param_name)
            
            elif current_section == 'shell':
                # Extract shell command
                if stripped_line.startswith('"') and stripped_line.endswith('"'):
                    metadata['shell_command'] = stripped_line.strip('"')
        
        logger.debug(f"Extracted metadata for rule '{rule_name}': {metadata}")
        return metadata
        
    except Exception as e:
        logger.error(f"Error extracting metadata for rule '{rule_name}': {e}")
        return None


def _get_available_rules(snakefile_path: str) -> List[str]:
    """Extract list of available rule names from a Snakefile."""
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Find all rule definitions
        rule_pattern = r'^rule\s+(\w+):'
        matches = re.findall(rule_pattern, content, re.MULTILINE)
        return matches
        
    except Exception as e:
        logger.error(f"Error getting available rules: {e}")
        return []


def _apply_replacements(content: str, parameters: Optional[Dict[str, Any]], 
                       file_mappings: Optional[Dict[str, str]]) -> str:
    """
    Apply parameter and file mapping replacements to extracted rule content.
    
    Args:
        content: The Snakefile rule content to process
        parameters: Dict of parameter name/value pairs to replace
        file_mappings: Dict mapping file placeholders to absolute paths
                      e.g., {"input.txt": "/absolute/path/to/expression.csv"}
                      These replace placeholders like {input.txt} with actual file paths
    
    Returns:
        str: Content with all placeholders replaced
        
    Raises:
        SnakefileGenerationError: If placeholder replacement fails
    """
    try:
        # Apply parameter replacements
        if parameters:
            for param_name, param_value in parameters.items():
                # Replace parameter placeholders like {Top Genes} or {Sample Size}
                # Handle various naming conventions - covers all cases without duplication
                placeholder_patterns = [
                    f"{{{param_name}}}",  # Exact match
                    f"{{{param_name.replace('_', ' ')}}}",  # Underscore to space
                    f"{{{param_name.replace(' ', '_')}}}",  # Space to underscore
                ]
                
                for pattern in placeholder_patterns:
                    if pattern in content:
                        content = content.replace(pattern, str(param_value))
                        logger.debug(f"Replaced parameter placeholder: {pattern} -> {param_value}")
        
        # Apply file mapping replacements
        if file_mappings:
            for placeholder, actual_path in file_mappings.items():
                # Replace file placeholders in input/output sections
                placeholder_patterns = [
                    f"{{{placeholder}}}",  # Standard placeholder format
                    f"/{{{placeholder}}}",  # With leading slash
                    f'"{{{placeholder}}}"',  # In quotes
                    f"'{{{placeholder}}}'",  # In single quotes
                ]
                
                for pattern in placeholder_patterns:
                    if pattern in content:
                        # Preserve quotes if they exist
                        if pattern.startswith('"') and pattern.endswith('"'):
                            replacement = f'"{actual_path}"'
                        elif pattern.startswith("'") and pattern.endswith("'"):
                            replacement = f"'{actual_path}'"
                        else:
                            replacement = actual_path
                        
                        content = content.replace(pattern, replacement)
                        logger.debug(f"Replaced file placeholder: {pattern} -> {replacement}")
        
        # Check for any remaining placeholders that weren't replaced
        remaining_placeholders = re.findall(r'\{[^}]+\}', content)
        if remaining_placeholders:
            logger.warning(f"Unreplaced placeholders found: {remaining_placeholders}")
            logger.debug(f"Available parameters: {list(parameters.keys()) if parameters else []}")
            logger.debug(f"Available file mappings: {list(file_mappings.keys()) if file_mappings else []}")
        
        return content
        
    except Exception as e:
        logger.error(f"Error applying replacements: {e}")
        # Re-raise the exception to prevent silent failures
        raise SnakefileGenerationError(
            plugin_name="unknown",
            script_name="unknown", 
            error_details=f"Placeholder replacement failed: {str(e)}"
        ) from e


def _extract_filename_from_pattern(file_pattern: str) -> Optional[str]:
    """Extract filename from Snakemake file pattern."""
    try:
        # Handle patterns like "user/{user_name}/data/{filename}"
        if "/" in file_pattern:
            # Get the last part of the path
            filename = file_pattern.split("/")[-1]
        else:
            filename = file_pattern
        
        # Remove curly braces if present
        filename = filename.strip("{}")
        
        # Return None for template variables
        if filename in ['user_name', 'workflow_id', 'algorithm_id']:
            return None
            
        return filename
        
    except Exception:
        return None