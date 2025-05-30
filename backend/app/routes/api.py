from fastapi import APIRouter

from app.routes.endpoints import auth, workflow, files, admin, plugin, task, datatable

api_router = APIRouter()

api_router.include_router(auth.router, prefix="/auth", tags=["auth"])
api_router.include_router(workflow.router, prefix="/workflow", tags=["workflow"])
api_router.include_router(files.router, prefix="/files", tags=["files"])
api_router.include_router(admin.router, prefix="/admin", tags=["admin"])
api_router.include_router(plugin.router, prefix="/plugin", tags=["plugin"])
api_router.include_router(task.router, prefix="/task", tags=["task"])
api_router.include_router(datatable.router, prefix="/datatable", tags=["datatable"])