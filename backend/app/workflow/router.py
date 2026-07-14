from fastapi import APIRouter, Depends
from typing import Any
from sqlalchemy.orm import Session
import logging

from app.workflow.schemas import (
    WorkflowDelete, WorkflowCreate, WorkflowResult, WorkflowFind,
    WorkflowNodeFileCreate, WorkflowNodeFileDelete, WorkflowNodeFileRead,
    WorkflowVisualizationRequest
)
from app.auth import deps as dep
from app import models
from app.workflow import service

# --- Test-patch compatibility re-exports ---------------------------------
# The business logic now lives in ``app.workflow.service``; the router only
# wires HTTP concerns to it. Existing tests (test_workflow_api.py,
# test_characterization_workflow_compile.py) patch these symbols on the router
# namespace (e.g. ``app.workflow.router.get_plugin_path``). Keeping them
# importable here preserves those patch targets without changing behavior.
# noqa: F401 imports below are intentional re-exports for test compatibility.
from app.plugin.utils import get_plugin_path  # noqa: F401
from app.worker.utils import get_task_info  # noqa: F401
from app.workflow.compiler.snakefile import change_snakefile_parameter  # noqa: F401
from app.workflow.utils import extract_rule_block  # noqa: F401
from app.workflow import crud as crud_workflow  # noqa: F401
from app.plugin import crud as crud_plugin  # noqa: F401
from app.worker.tasks import process_data_task  # noqa: F401

router = APIRouter()
logger = logging.getLogger(__name__)


# compile workflow
@router.post("/compile")
def compileWorkflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.compile_workflow(db=db, user=current_user, workflow=workflow)


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
    return service.visualize_data(db=db, user=current_user, workflow_request=workflow_request)


# get visualization result
@router.post("/visualize/result")
def getVisualizationResult(
    WorkflowResult: WorkflowResult,
    current_user: models.User = Depends(dep.get_current_active_user)
):
    return service.get_visualization_result(workflow_result=WorkflowResult)


# save workflow data
@router.post("/save")
def update_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.save_workflow(db=db, user=current_user, workflow=workflow)


# User workflow delete
@router.post("/delete")
def delete_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowDelete,
    ) -> Any:
    return service.delete_workflow(db=db, user=current_user, workflow_id=workflowInfo.id)


# response workflow data
@router.get("/me")
def get_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    return service.list_workflows(db=db, user=current_user)


# find user workflow
@router.post("/find", response_model=WorkflowCreate)
def find_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowFind,
    ) -> Any:
    return service.find_workflow(db=db, user=current_user, workflow_id=workflowInfo.id)


# response workflow results
@router.post("/results")
def getResults(
    WorkflowResult: WorkflowResult,
    current_user: models.User = Depends(dep.get_current_active_user)
):
    return service.get_results(user=current_user, workflow_result=WorkflowResult)


# response workflow result
@router.post("/result")
def checkResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    return service.check_result(user=current_user, workflow_result=WorkflowResult)


# response workflow visualization result
@router.post("/visualization/result")
def checkVisualizationResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    return service.check_visualization_result(user=current_user, workflow_result=WorkflowResult)


# save workflow node modal data
@router.post("/node/save")
def saveNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileCreate,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.save_node_data(db=db, user=current_user, node_file_info=workflowNodeFileInfo)


# read workflow node modal data
@router.post("/node/read")
def readNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileRead,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.read_node_data(db=db, user=current_user, node_file_info=workflowNodeFileInfo)


# delete workflow node modal data
@router.post("/node/delete")
def deleteNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileDelete,
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    return service.delete_node_data(db=db, user=current_user, node_file_info=workflowNodeFileInfo)
