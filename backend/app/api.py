from fastapi import APIRouter

from app.auth.router import router as auth_router
from app.workflow.router import router as workflow_router
from app.file.router import router as files_router
from app.admin.router import router as admin_router
from app.plugin.router import router as plugin_router
from app.task.router import router as task_router
from app.datatable.router import router as datatable_router
from app.version.router import router as version_router

api_router = APIRouter()

api_router.include_router(auth_router, prefix="/auth", tags=["auth"])
api_router.include_router(workflow_router, prefix="/workflow", tags=["workflow"])
api_router.include_router(files_router, prefix="/files", tags=["files"])
api_router.include_router(admin_router, prefix="/admin", tags=["admin"])
api_router.include_router(plugin_router, prefix="/plugin", tags=["plugin"])
api_router.include_router(task_router, prefix="/task", tags=["task"])
api_router.include_router(datatable_router, prefix="/datatable", tags=["datatable"])
api_router.include_router(version_router, prefix="/version", tags=["version"])