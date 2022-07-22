from fastapi import APIRouter

from app.routes.endpoints import auth, workflow

api_router = APIRouter()

api_router.include_router(auth.router, prefix="/auth", tags=["auth"])
api_router.include_router(workflow.router, prefix="/workflow", tags=["workflow"])