from fastapi import FastAPI
from sqlalchemy import engine
from starlette.middleware.cors import CORSMiddleware

from app.routes.api import api_router
from app.common.config import settings
from app.common.utils.celery_utils import create_celery
from app.database import models
from app.database.conn import engine
from celery.signals import worker_shutting_down

# Celery 이벤트 핸들러를 정의합니다
def on_worker_shut_down(sender=None, conf=None, **kwargs):
    print("Worker is shutting down. Trying to reconnect...")
    try:
        with celery.connection() as connection:
            connection.ensure_connection(max_retries=3)
            print("재연결에 성공했습니다.")
    except Exception as e:
        print(f"재연결 시도 중 오류 발생: {e}")

models.Base.metadata.create_all(bind=engine)
global_engine = engine

app = FastAPI(
    title=settings.PROJECT_NAME
)

app.add_middleware(
    CORSMiddleware,
    # allow_origins=[str(origin) for origin in settings.BACKEND_CORS_ORIGINS],
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.celery_app = create_celery()
# Celery 이벤트 핸들러를 Celery 애플리케이션 인스턴스에 연결합니다
worker_shutting_down.connect(on_worker_shut_down, sender=app.celery_app)

app.include_router(api_router, prefix=settings.ROUTES_STR)

celery = app.celery_app

def get_celery_app():
    return celery