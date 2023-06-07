from fastapi import FastAPI
from sqlalchemy import engine
from starlette.middleware.cors import CORSMiddleware

from app.routes.api import api_router
from app.common.config import settings
from app.common.celery_utils import create_celery
from app.database import models
from app.database.conn import engine


models.Base.metadata.create_all(bind=engine)

global_engine = engine

app = FastAPI(
    title=settings.PROJECT_NAME
)

origins = [
    "http://localhost.tiangolo.com",
    "https://localhost.tiangolo.com",
    "http://localhost",
    "http://localhost:8080",
]

# if settings.BACKEND_CORS_ORIGINS:

app.add_middleware(
    CORSMiddleware,
    # allow_origins=[str(origin) for origin in settings.BACKEND_CORS_ORIGINS],
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.celery_app = create_celery()
app.include_router(api_router, prefix=settings.ROUTES_STR)

celery = app.celery_app