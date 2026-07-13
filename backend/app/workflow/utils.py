"""Workflow 유틸리티 (PR-10 부분 분할).

compile 파이프라인의 Snakefile rule 추출/생성 로직은 ``compiler/visualization.py``로,
공유 예외 계층은 ``workflow/errors.py``로 이동했다. 이 모듈은 데이터 로딩/그래프
탐색/파라미터 생성/파일 경로 검증 등 잔류 헬퍼를 담고, 이동한 심볼들을 아래에서
그대로 re-export 하여 기존 import 경로(``app.workflow.utils.extract_rule_block``,
``app.workflow.utils.FileValidationError`` 등)를 유지한다.
"""
import os
import pandas as pd
import scanpy as sc
import numpy as np
import logging
import json
import gc
from typing import Dict, List, Optional, Tuple, Any

# Configure logging
logger = logging.getLogger(__name__)

# --- Facade re-exports (PR-10 split) -------------------------------------
# 예외 계층은 workflow/errors.py로, Snakefile 처리 함수는 compiler/visualization.py로
# 이동했다. 기존 참조 경로 호환을 위해 여기서 re-export 한다.
from app.workflow.errors import (  # noqa: F401
    WorkflowUtilsError,
    FileValidationError,
    SnakefileProcessingError,
    AlgorithmResolutionError,
)
from app.workflow.compiler.visualization import (  # noqa: F401
    extract_rule_block,
    generate_visualization_snakefile,
    validate_visualization_rule,
    get_rule_metadata,
    _get_available_rules,
    _apply_replacements,
    _extract_filename_from_pattern,
)
# -------------------------------------------------------------------------

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
