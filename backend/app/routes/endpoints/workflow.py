from fastapi import APIRouter, Depends, HTTPException
from typing import Any
from sqlalchemy.orm import Session
import os
import shutil
import json
import logging
import time
from pathlib import Path
from fastapi.responses import FileResponse, JSONResponse

from app.common.utils.plugin_utils import verify_dependencies, get_plugin_path
from app.common.utils.celery_utils import get_task_info
from app.common.utils.snakemake_utils import change_snakefile_parameter
from app.common.utils.error_utils import (
    PluginNotFoundError, ScriptNotFoundError, FileNotFoundError as CellCraftFileNotFoundError,
    ValidationError, WorkflowError, SnakefileGenerationError, TaskSubmissionError,
    log_error, create_error_response
)
from app.common.utils.cache_utils import (
    generate_cache_key, check_cache_with_expiry, create_symbolic_link,
    save_result_to_cache, maybe_cleanup_cache, update_cache_link_location,
    remove_cache_by_visualization_path
)
from app.common.utils.workflow_utils import (
    extract_rule_block, extract_all_algorithms, extract_algorithm_data,
    extract_visualization_data, extract_target_data, generate_user_input,
    generate_plugin_params, generate_visualization_params,
    resolve_algorithm_path_from_files, validate_file_paths,
    find_connected_visualization_nodes, extract_file_sources
)
from app.common.utils.log_archive_utils import cleanup_task_results
from app.database.crud import crud_workflow, crud_plugin
from app.database.schemas.workflow import (
    WorkflowDelete, WorkflowCreate, WorkflowUpdate, WorkflowResult, WorkflowFind, 
    WorkflowNodeFileCreate, WorkflowNodeFileDelete, WorkflowNodeFileRead, 
    WorkflowVisualizationRequest
)
from app.routes import dep
from app.database import models
from app.routes.celery_tasks import process_data_task

router = APIRouter()
logger = logging.getLogger(__name__)

#compile workflow
@router.post("/compile")
def compileWorkflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        user_path = f"./user/{current_user.username}/"
        task_ids = []  # 여러 개의 태스크 ID를 저장할 리스트

        if user_workflow:
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info)
            algorithms = extract_all_algorithms(user_workflow.workflow_info['drawflow']['Home']['data'])

            if not algorithms:
                raise HTTPException(status_code=400, detail="No algorithm nodes found in workflow.")

            for algorithm in algorithms:
                selected_plugin_name = algorithm['selectedPlugin']['name']
                selected_plugin_source = algorithm['selectedPlugin'].get('source')  # Get source from frontend
                user_workflow_task_path = f"{user_path}workflow_{workflow.id}/algorithm_{algorithm['id']}"

                # Get plugin path based on its source (official/local)
                try:
                    if selected_plugin_source:
                        # Use source information from frontend if available
                        plugin_path, actual_source = get_plugin_path(selected_plugin_name, selected_plugin_source)
                        print(f"Using plugin source from frontend: {selected_plugin_source} -> {actual_source}")
                    else:
                        # Fallback to auto-detection via database lookup
                        plugin_path, actual_source = get_plugin_path(selected_plugin_name)
                        print(f"Auto-detected plugin source: {actual_source}")
                    
                    plugin_dependency_path = os.path.join(plugin_path, "dependency")
                except Exception as e:
                    raise HTTPException(status_code=404, detail=f"Plugin '{selected_plugin_name}' not found: {str(e)}")

                # 입력 데이터 및 파라미터 추출 (파일 소스에 따라 전체 경로 구성)
                file_sources = extract_file_sources(algorithm)
                user_input = generate_user_input(
                    algorithm['selectedPluginInputOutput'],
                    username=current_user.username,
                    file_sources=file_sources,
                )
                plugin_params = generate_plugin_params(algorithm['selectedPluginRules'])
                target_list = extract_target_data(algorithm['selectedPluginInputOutput'], user_workflow_task_path)

                additional_data = {
                    "user_name": current_user.username,
                    "workflow_id": str(workflow.id),
                    "algorithm_id": str(algorithm['id']),
                    "plugin_name": selected_plugin_name,
                }
                user_input.update(additional_data)

                # 작업 폴더 생성 또는 기존 results 정리
                if not os.path.exists(user_workflow_task_path):
                    os.makedirs(user_workflow_task_path)
                else:
                    # 재실행 시 기존 results 폴더 정리 (Snakemake 스킵 방지)
                    cleanup_result = cleanup_task_results(Path(user_workflow_task_path), preserve_folder=True)
                    if cleanup_result["success"]:
                        files_count = len(cleanup_result.get('files_removed', [])) + len(cleanup_result.get('symlinks_removed', []))
                        if files_count > 0:
                            logger.info(f"Cleaned up {files_count} previous result files for algorithm_{algorithm['id']}")

                    # 연결된 Visualization 노드 폴더도 정리
                    workflow_data = user_workflow.workflow_info['drawflow']['Home']['data']
                    connected_vis_ids = find_connected_visualization_nodes(workflow_data, algorithm['id'])
                    for vis_id in connected_vis_ids:
                        vis_path = Path(f"{user_path}workflow_{workflow.id}/visualization_{vis_id}")
                        if vis_path.exists():
                            # 캐시 정리 먼저 수행 (symlink 삭제 전에 메타데이터 참조 필요)
                            vis_relative_path = f"workflow_{workflow.id}/visualization_{vis_id}/results"
                            cache_cleanup = remove_cache_by_visualization_path(user_path, vis_relative_path)
                            if cache_cleanup.get("cache_files_removed"):
                                logger.info(f"Cleaned up {len(cache_cleanup['cache_files_removed'])} cache files for visualization_{vis_id}")

                            # 기존 결과 폴더 정리 (symlinks 및 파일)
                            vis_cleanup = cleanup_task_results(vis_path, preserve_folder=True)
                            if vis_cleanup["success"]:
                                vis_files_count = len(vis_cleanup.get('files_removed', [])) + len(vis_cleanup.get('symlinks_removed', []))
                                if vis_files_count > 0:
                                    logger.info(f"Cleaned up {vis_files_count} previous files for visualization_{vis_id}")

                # Snakefile 생성
                plugin_snakefile_path = os.path.join(plugin_path, "Snakefile")
                user_snakefile_path = change_snakefile_parameter(plugin_snakefile_path, user_workflow_task_path + "/Snakefile", user_input, plugin_params)

                # Celery 작업 실행
                process_task = process_data_task.apply_async(
                    args=[current_user.username, user_snakefile_path, selected_plugin_name, target_list],
                    kwargs={'user_id': current_user.id, 'workflow_id': workflow.id, 'algorithm_id': algorithm['id'], 'plugin_name': selected_plugin_name, 'task_type': 'compile'},
                    ignore_result=False
                )

                # 실행된 태스크 ID 저장
                task_ids.append(process_task.id)

            # Extract algorithm IDs for frontend state management
            algorithm_ids = [algorithm['id'] for algorithm in algorithms]
            task_algorithm_mapping = dict(zip(task_ids, algorithm_ids))
            
            return {
                "message": "Multiple tasks added to queue",
                "task_ids": task_ids,
                "algorithm_ids": algorithm_ids,
                "task_algorithm_mapping": task_algorithm_mapping,
                "results": [get_task_info(task_id) for task_id in task_ids]
            }

    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
        )

# visualize compile
@router.post("/visualization")
def visualizeData(
    *,
    db: Session = Depends(dep.get_db),
    workflow_request: WorkflowVisualizationRequest, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    """
    Create visualization from workflow data using the new self-contained visualization system.
    
    This endpoint processes visualization requests using:
    - selectedVisualizationPlugin to get plugin information
    - selectedScript.name to identify the specific visualization rule
    - resolve_algorithm_path_from_files() for file resolution
    - generate_visualization_snakefile() for Snakefile generation
    
    Raises:
        ValidationError: If required parameters are missing or invalid
        WorkflowError: If workflow is not found or inaccessible
        PluginNotFoundError: If specified plugin is not available
        ScriptNotFoundError: If visualization script is not found
        FileNotFoundError: If required files are missing
        SnakefileGenerationError: If Snakefile generation fails
        TaskSubmissionError: If task submission fails
    """
    logger.info(f"Processing visualization request for workflow {workflow_request.id}")
    
    try:
        # Get user workflow
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow_request.id)
        if not user_workflow:
            raise WorkflowError(
                workflow_id=str(workflow_request.id),
                message=f"Workflow with ID {workflow_request.id} not found or not accessible by user {current_user.username}",
                error_type="not_found"
            )
        
        # Update workflow if additional data provided
        if workflow_request.workflow_info:
            crud_workflow.update_workflow(
                db, current_user.id, workflow_request.id, 
                workflow_request.title, workflow_request.thumbnail, 
                workflow_request.workflow_info
            )
        
        # Extract plugin information from selectedVisualizationPlugin
        selected_plugin = workflow_request.selectedVisualizationPlugin
        if not selected_plugin:
            raise ValidationError(
                field_name="selectedVisualizationPlugin",
                message="Plugin selection is required for visualization",
                required_format="{name: string, source?: string}"
            )
        
        selected_plugin_name = selected_plugin.get('name')
        selected_plugin_source = selected_plugin.get('source')
        
        if not selected_plugin_name:
            raise ValidationError(
                field_name="selectedVisualizationPlugin.name",
                message="Plugin name is required",
                required_format="string"
            )
        
        # Extract visualization script name from selectedScript
        selected_script = workflow_request.selectedScript
        if not selected_script:
            raise ValidationError(
                field_name="selectedScript",
                message="Script selection is required for visualization",
                required_format="{name: string, parameters?: array}"
            )
        
        selected_script_name = selected_script.get('name')
        if not selected_script_name:
            raise ValidationError(
                field_name="selectedScript.name",
                message="Script name is required",
                required_format="string"
            )
        
        logger.info(f"Processing visualization: plugin={selected_plugin_name}, script={selected_script_name}")
        
        # Get plugin path based on source
        try:
            if selected_plugin_source:
                plugin_path, actual_source = get_plugin_path(selected_plugin_name, selected_plugin_source)
                logger.info(f"Using plugin source from request: {selected_plugin_source} -> {actual_source}")
            else:
                plugin_path, actual_source = get_plugin_path(selected_plugin_name)
                logger.info(f"Auto-detected plugin source: {actual_source}")
        except Exception as e:
            log_error(e, {"plugin_name": selected_plugin_name, "source": selected_plugin_source})
            # Try to get available plugins for better error message
            try:
                available_plugins = crud_plugin.get_available_plugin_names(db)
            except Exception:
                available_plugins = None
            raise PluginNotFoundError(selected_plugin_name, available_plugins)
        
        # Set up paths
        user_path = f"./user/{current_user.username}/"
        workflow_info = workflow_request.workflow_info or user_workflow.workflow_info
        
        # Resolve algorithm path using new method
        algorithm_id = None
        if workflow_request.algorithm_id:
            # Use provided algorithm_id
            algorithm_id = workflow_request.algorithm_id
            logger.info(f"Using provided algorithm_id: {algorithm_id}")
        elif workflow_request.availableFiles:
            # Resolve from workflow connections and available files
            try:
                algorithm_id = resolve_algorithm_path_from_files(
                    workflow_info['drawflow']['Home']['data'],
                    workflow_request.current_node_id,
                    workflow_request.availableFiles
                )
                logger.info(f"Resolved algorithm_id from files: {algorithm_id}")
            except Exception as e:
                log_error(e, {
                    "workflow_id": workflow_request.id,
                    "current_node_id": workflow_request.current_node_id
                })
                raise ValidationError(
                    field_name="algorithm_id",
                    message=f"Could not resolve algorithm path: {str(e)}",
                    required_format="Valid algorithm node connection or explicit algorithm_id"
                )
        
        if not algorithm_id:
            raise ValidationError(
                field_name="algorithm_id",
                message="Cannot determine algorithm_id for file resolution",
                required_format="Either provide algorithm_id or ensure proper workflow connections"
            )
        
        # Set up task paths for visualization
        user_task_path = f"{user_path}workflow_{workflow_request.id}/visualization_{workflow_request.current_node_id}"
        user_workflow_visualization_result_directory = f"{user_task_path}/results"
        
        # Process visualization parameters
        visualization_params = {}
        file_mappings = {}
        input_filenames = []
        output_filename = None
        
        for param in workflow_request.selectedVisualizationParams:
            param_type = param.get("type")
            param_name = param.get("name")
            param_value = param.get("defaultValue")
            
            if param_type in ["inputFile", "optionalInputFile"]:
                # Input file parameter - get the actual selected file
                # selectedFile contains the user's selection, defaultValue contains the placeholder name
                selected_file = param.get("selectedFile")
                placeholder_name = param.get("defaultValue", param_name)  # This is the placeholder like "input.txt"
                
                if selected_file:
                    input_filenames.append(selected_file)
                    # Map placeholder name to actual selected file
                    file_mappings[placeholder_name] = selected_file
                elif param_value:
                    # Fallback if selectedFile is not provided but defaultValue is
                    input_filenames.append(param_value)
                    file_mappings[placeholder_name] = param_value
            elif param_type == "outputFile":
                # Output file parameter
                output_filename = param_value
            else:
                # Regular parameter
                visualization_params[param_name] = param_value
        
        logger.debug(f"Visualization params: {visualization_params}")
        logger.debug(f"File mappings: {file_mappings}")
        logger.debug(f"Input filenames: {input_filenames}")
        logger.debug(f"Output filename: {output_filename}")
        
        # Validate file paths if we have input files
        validated_file_paths = []
        if input_filenames:
            try:
                validated_file_paths = validate_file_paths(
                    user_path, workflow_request.id, algorithm_id, input_filenames
                )
                logger.info(f"Validated {len(validated_file_paths)} file paths")
                
                # Update file mappings with just the filenames (not full paths)
                # The Snakefile template already has the correct base path structure
                for i, filename in enumerate(input_filenames):
                    if i < len(validated_file_paths):
                        # Find the placeholder key that maps to this filename
                        for placeholder_key, mapped_filename in file_mappings.items():
                            if mapped_filename == filename:
                                # Use only the basename to avoid path duplication
                                # The template already has the correct directory structure
                                file_mappings[placeholder_key] = os.path.basename(validated_file_paths[i])
                                break
                                
            except FileNotFoundError as e:
                log_error(e, {
                    "user_path": user_path,
                    "workflow_id": workflow_request.id,
                    "algorithm_id": algorithm_id,
                    "input_filenames": input_filenames
                })
                raise CellCraftFileNotFoundError(
                    file_path=f"algorithm_{algorithm_id}/results",
                    file_type="input files"
                )
            except ValueError as e:
                log_error(e, {"input_filenames": input_filenames})
                raise ValidationError(
                    field_name="input_files",
                    message=f"File validation failed: {str(e)}",
                    required_format="Valid file paths in algorithm results directory"
                )
            except Exception as e:
                log_error(e)
                raise ValidationError(
                    field_name="input_files",
                    message=f"Unexpected error during file validation: {str(e)}"
                )
        
        # Generate cache key for per-node caching (includes workflow_id and node_id for isolation)
        cache_key = generate_cache_key(
            selected_plugin_name,
            selected_script_name,
            visualization_params,
            input_filenames,
            workflow_id=workflow_request.id,
            node_id=workflow_request.current_node_id
        )
        
        # Run probabilistic cache cleanup (10% chance)
        maybe_cleanup_cache(user_path)
        
        # Generate result filename based on parameters and input files (original logic)
        result_name_parts = []
        for key, value in visualization_params.items():
            # Normalize parameter key by replacing spaces with underscores
            # This matches the normalize_param_name behavior in Snakefile generation
            normalized_key = key.replace(" ", "_")
            result_name_parts.append(f"{normalized_key}_{value}")
        for filename in input_filenames:
            result_name_parts.append(os.path.basename(filename))
        
        result_name = "_".join(result_name_parts)
        
        # Set up result file path with parameter-based naming
        if output_filename:
            if result_name_parts:  # If we have parameters, prepend them
                result_filename = f"{result_name}_{output_filename}"
            else:
                result_filename = output_filename
        else:
            # Default output filename if not specified
            if result_name_parts:
                result_filename = f"{result_name}_visualization_result.json"
            else:
                result_filename = "visualization_result.json"
            
        user_workflow_visualization_result_path = os.path.join(
            user_workflow_visualization_result_directory, 
            result_filename
        )
        
        # Check centralized cache with expiry
        cache_hit, cache_file_path = check_cache_with_expiry(user_path, cache_key)
        
        if cache_hit:
            logger.info(f"Cache hit for key {cache_key}: {cache_file_path}")
            
            # Create task directory if it doesn't exist
            os.makedirs(user_task_path, exist_ok=True)
            os.makedirs(user_workflow_visualization_result_directory, exist_ok=True)
            
            # Create symbolic link to cached result
            if create_symbolic_link(cache_file_path, user_workflow_visualization_result_path):
                # Update cache metadata with new link location
                relative_link_path = os.path.relpath(
                    user_workflow_visualization_result_path, 
                    user_path
                )
                update_cache_link_location(user_path, cache_key, relative_link_path)
                
                return {
                    "message": "Visualization result from cache",
                    "result_path": result_filename,
                    "cached": True
                }
            else:
                logger.warning(f"Failed to create symbolic link for cache key {cache_key}")
                # Continue with normal execution if link creation fails
        
        # Create task directory if it doesn't exist
        os.makedirs(user_task_path, exist_ok=True)
        
        # Load plugin metadata to get rules data
        try:
            plugin_info = crud_plugin.get_plugin_by_name(db, selected_plugin_name)
            if not plugin_info:
                raise PluginNotFoundError(selected_plugin_name)
            
            # Get plugin metadata
            metadata_path = os.path.join(plugin_path, "metadata.json")
            with open(metadata_path, 'r') as f:
                plugin_metadata = json.load(f)
            
            rules_data = plugin_metadata.get('rules', {})
            if not rules_data:
                raise ValueError(f"No rules found in plugin {selected_plugin_name}")
            
            # Filter to get only the selected visualization script
            selected_rule = None
            for rule_id, rule in rules_data.items():
                if rule.get('name') == selected_script_name:
                    selected_rule = {rule_id: rule}
                    break
            
            if not selected_rule:
                raise ScriptNotFoundError(selected_script_name, selected_plugin_name)
            
            # Generate final Snakefile directly from template
            plugin_snakefile_path = os.path.join(plugin_path, "Snakefile")
            user_snakefile_path = os.path.join(user_task_path, "Snakefile")
            
            logger.info(f"Using template Snakefile: {plugin_snakefile_path}")
            logger.info(f"Generating final Snakefile: {user_snakefile_path}")
            
        except FileNotFoundError as e:
            log_error(e, {"plugin_path": plugin_path, "metadata_path": metadata_path})
            raise PluginNotFoundError(selected_plugin_name)
        except json.JSONDecodeError as e:
            log_error(e, {"metadata_path": metadata_path})
            raise ValidationError(
                field_name="plugin_metadata",
                message=f"Invalid plugin metadata format: {str(e)}"
            )
        except Exception as e:
            log_error(e, {
                "plugin_path": plugin_path,
                "script_name": selected_script_name,
                "parameters": visualization_params
            })
            raise SnakefileGenerationError(
                plugin_name=selected_plugin_name,
                script_name=selected_script_name,
                error_details=str(e)
            )
        
        # Additional context data for task execution
        additional_data = {
            "user_name": current_user.username,
            "workflow_id": str(workflow_request.id),
            "algorithm_id": str(algorithm_id),
            "visualization_id": str(workflow_request.current_node_id),
            "plugin_name": selected_plugin_name,
            "visualization_name": selected_script_name,
            "result_filename": result_filename  # Add the dynamic filename
        }
        
        # file_mappings contains placeholder-to-path mappings for replacement in the Snakefile
        # e.g., {"input.txt": "/absolute/path/to/expression.csv", "trajectory.txt": "/absolute/path/to/trajectory.csv"}
        # These will replace placeholders like {input.txt} with actual file paths
        
        # Extract only the selected rule and apply replacements
        try:
            # Merge system placeholders with visualization parameters
            all_params = {**visualization_params, **additional_data}
            
            # Use extract_rule_block to get only the selected visualization rule
            extracted_content, final_snakefile_path = extract_rule_block(
                snakefile_path=plugin_snakefile_path,
                visualization_title=selected_script_name,
                result_path=user_snakefile_path,
                parameters=all_params,  # Now includes both user params and system placeholders
                file_mappings=file_mappings
            )
            logger.info(f"Successfully extracted and generated rule-specific Snakefile: {final_snakefile_path}")
        except Exception as e:
            logger.error(f"Failed to extract rule from Snakefile: {e}")
            raise SnakefileGenerationError(
                plugin_name=selected_plugin_name,
                script_name=selected_script_name,
                error_details=f"Rule extraction failed: {str(e)}"
            )
        
        # Set target list for Celery task
        target_list = [user_workflow_visualization_result_path]
        
        # Execute Celery task
        try:
            process_task = process_data_task.apply_async(
                args=[current_user.username, final_snakefile_path, selected_plugin_name, target_list],
                kwargs={
                    'user_id': current_user.id,
                    'workflow_id': workflow_request.id,
                    'algorithm_id': workflow_request.current_node_id,  # Use visualization node ID for correct log path
                    'plugin_name': selected_plugin_name,
                    'task_type': 'visualization',
                    'cache_key': cache_key,
                    'cache_info': {
                        'user_path': user_path,
                        'result_file_path': user_workflow_visualization_result_path,
                        'plugin_name': selected_plugin_name,
                        'script_name': selected_script_name,
                        'linked_location': os.path.relpath(user_workflow_visualization_result_path, user_path)
                    }
                },
                ignore_result=False
            )
            
            task_id = process_task.id
            task_info = get_task_info(task_id)
            
            logger.info(f"Visualization task submitted: {task_id}")
            
            return {
                "message": "Visualization task added to queue",
                "task_id": task_id,
                "result": task_info,
                "result_path": result_filename,
                "cached": False
            }
            
        except Exception as e:
            log_error(e, {
                "user": current_user.username,
                "plugin_name": selected_plugin_name,
                "snakefile_path": final_snakefile_path
            })
            raise TaskSubmissionError(
                task_type="visualization",
                error_details=str(e)
            )
        
    except (PluginNotFoundError, ScriptNotFoundError, CellCraftFileNotFoundError, 
            ValidationError, WorkflowError, SnakefileGenerationError, TaskSubmissionError):
        # Re-raise our custom HTTP exceptions as-is
        raise
    except HTTPException:
        # Re-raise standard HTTP exceptions as-is
        raise
    except Exception as e:
        log_error(e, {"workflow_id": workflow_request.id, "user": current_user.username})
        # Create a structured error response for unexpected errors
        return create_error_response(
            error=e,
            default_message="An unexpected error occurred during visualization processing"
        )

# get visualization result
@router.post("/visualize/result")
def getVisualizationResult(
    WorkflowResult: WorkflowResult, 
    current_user: models.User = Depends(dep.get_current_active_user)
):
    PATH_VISUALIZAE_RESULT = WorkflowResult.filename
    print("PATH_VISUALIZAE_RESULT:", PATH_VISUALIZAE_RESULT)

    try:
        with open(PATH_VISUALIZAE_RESULT, "r") as f:
            plotly_data = json.load(f)
        return JSONResponse(content=plotly_data)
        
    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
                )

#save workflow data
@router.post("/save")
def update_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
    if user_workflow:
        # workflow 수정
        return crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info)
    else :
        # workflow 생성
        return crud_workflow.create_workflow(db, workflow.title, workflow.thumbnail, workflow.workflow_info, current_user.id)

#User workflow delete
@router.post("/delete")
def delete_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowDelete,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        # 실제 폴더 경로 구성 (언더스코어 추가)
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowInfo.id}"
        
        # 보안: Path traversal 방지
        real_user_path = os.path.realpath(user_path)
        real_workflow_path = os.path.realpath(workflow_path)
        
        # 워크플로우 경로가 사용자 폴더 내에 있는지 확인
        if not real_workflow_path.startswith(real_user_path):
            raise HTTPException(
                status_code=403,
                detail="Access denied: Invalid workflow path"
            )
        
        # 폴더 존재 여부 확인 및 삭제
        if os.path.exists(workflow_path) and os.path.isdir(workflow_path):
            try:
                shutil.rmtree(workflow_path)  # 폴더 및 하위 파일 모두 삭제
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to delete workflow folder: {str(e)}"
                )
        
        # 폴더 삭제 성공 후 DB에서 워크플로우 정보 삭제
        delete_workflow = crud_workflow.delete_user_workflow(db, current_user.id, workflowInfo.id)
        return delete_workflow
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#response workflow data
@router.get("/me")
def get_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    user_workflows = crud_workflow.get_user_workflows(db, current_user.id)
    if user_workflows:
        res = []
        for item in user_workflows:
            workflow_res = {
                'id': item.id, 
                'title': item.title, 
                'thumbnail': item.thumbnail,
                'updated_at': item.updated_at, 
                'user_id': item.user_id
            }
            res.append(workflow_res)
        # print(res)
        return res
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#find user workflow
@router.post("/find", response_model=WorkflowCreate)
def find_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowFind,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        return {
            'title': user_workflow.title,
            'thumbnail': user_workflow.thumbnail,
            'workflow_info': user_workflow.workflow_info,
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#response workflow results
@router.post("/results")
def getResults(
    WorkflowResult: WorkflowResult, 
    current_user: models.User = Depends(dep.get_current_active_user)
):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'
    
    # 파일 리스트와 각 파일의 부가 정보 가져오기
    file_info_list = []
    if os.path.exists(PATH_COMPILE_RESULT):
        for file_name in os.listdir(PATH_COMPILE_RESULT):
            file_path = os.path.join(PATH_COMPILE_RESULT, file_name)
            if os.path.isfile(file_path):
                file_info = {
                    "name": file_name,
                    "size": os.path.getsize(file_path),  # 파일 크기
                    "modified_time": os.path.getmtime(file_path)  # 마지막 수정 시간
                }
                file_info_list.append(file_info)
    else:
        raise HTTPException(status_code=404, detail="Results directory not found.")
    
    return file_info_list

#response workflow result
@router.post("/result")
def checkResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'
    file_list = os.listdir(PATH_COMPILE_RESULT)
    # print(file_list)
    FILE_NAME = WorkflowResult.filename
    
    for item_file in file_list:
        if FILE_NAME == item_file:
            FILE_NAME = item_file
            break
    # print(FILE_NAME)
    FILE_PATH = os.path.join(PATH_COMPILE_RESULT, FILE_NAME)
    # print(FILE_PATH)

    return FileResponse(FILE_PATH)

#response workflow visualization result
@router.post("/visualization/result")
def checkVisualizationResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    # Check if this is a visualization node result (algorithm_id is actually a visualization node ID)
    # For visualization nodes, use visualization_{id} path instead of algorithm_{id}
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/visualization_{WorkflowResult.algorithm_id}/results'

    # If visualization path doesn't exist, try algorithm path for backward compatibility
    if not os.path.exists(PATH_COMPILE_RESULT):
        PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'

    file_list = os.listdir(PATH_COMPILE_RESULT)
    # print(file_list)
    FILE_NAME = WorkflowResult.filename

    for item_file in file_list:
        if FILE_NAME == item_file:
            FILE_NAME = item_file
            break
    # print(FILE_NAME)
    FILE_PATH = os.path.join(PATH_COMPILE_RESULT, FILE_NAME)
    # print(FILE_PATH)

    # Retry logic for handling Docker volume sync delays
    max_retries = 3
    retry_delay = 0.5  # 500ms between retries

    for attempt in range(max_retries):
        try:
            # Verify file exists and has content
            if os.path.exists(FILE_PATH) and os.path.getsize(FILE_PATH) > 0:
                with open(FILE_PATH, "r") as f:
                    plotly_data = json.load(f)

                # Log if retry was needed
                if attempt > 0:
                    logger.warning(f"Visualization file ready after {attempt} retry(ies): {FILE_PATH}")

                return JSONResponse(content=plotly_data)
            else:
                # File doesn't exist or is empty
                if attempt < max_retries - 1:
                    logger.info(f"Visualization file not ready, retry {attempt + 1}/{max_retries}: {FILE_PATH}")
                    time.sleep(retry_delay)
                else:
                    raise FileNotFoundError(f"File not found or empty: {FILE_PATH}")

        except (FileNotFoundError, json.JSONDecodeError) as e:
            if attempt < max_retries - 1:
                logger.info(f"Error reading visualization file, retry {attempt + 1}/{max_retries}: {str(e)}")
                time.sleep(retry_delay)
            else:
                logger.error(f"Failed to read visualization file after {max_retries} attempts: {FILE_PATH}")
                raise HTTPException(
                    status_code=400,
                    detail=str(e),
                )
        except Exception as e:
            # Unexpected error, don't retry
            logger.error(f"Unexpected error reading visualization file: {str(e)}")
            raise HTTPException(
                status_code=400,
                detail=str(e),
            )

#save workflow node modal data
@router.post("/node/save")
def saveNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # workflowNodeFileInfo.id로 사용자 폴더에 workflow 폴더 생성
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        # user 폴더에 workflow 폴더 존재하지 않으면 생성
        if not os.path.exists(workflow_path):
            os.makedirs(workflow_path)
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 이름 생성
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        result_file_path = f"{workflow_path}/{file_name}"
        # user 폴더에 파일 생성
        with open(result_file_path, "w") as f:
            # workflowNodeFileInfo.file_content가 List이면 json.dump으로 파일 생성
            if workflowNodeFileInfo.file_extension == "json":
                json.dump(workflowNodeFileInfo.file_content, f)
            if workflowNodeFileInfo.file_extension == "txt" or workflowNodeFileInfo.file_extension == "tsv" or workflowNodeFileInfo.file_extension == "csv":
                f.write(workflowNodeFileInfo.file_content)
        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_extension': workflowNodeFileInfo.file_extension,
            'file_path': result_file_path
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#read workflow node modal data
@router.post("/node/read")
def readNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileRead, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # user 폴더에 workflow 폴더 존재하지 않으면 에러 발생
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        print(workflow_path)
        if not os.path.exists(workflow_path):
            raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 읽어오기
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        with open(f"{workflow_path}/{file_name}", "r") as f:
            if workflowNodeFileInfo.file_extension == "json":
                file_content = json.load(f)
            else:
                file_content = f.read()

        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_content': file_content,
            'file_extension': workflowNodeFileInfo.file_extension
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#delete workflow node modal data
@router.post("/node/delete")
def deleteNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileDelete, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # user 폴더에 workflow 폴더 존재하지 않으면 에러 발생
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        if not os.path.exists(workflow_path):
            raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 삭제
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        os.remove(f"{workflow_path}/{file_name}")
        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_extension': workflowNodeFileInfo.file_extension
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
