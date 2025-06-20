from fastapi import FastAPI
from sqlalchemy import engine
from starlette.middleware.cors import CORSMiddleware
import signal
import sys

from app.routes.api import api_router
from app.common.config import settings
from app.common.utils.celery_utils import create_celery
from app.common.utils.docker_utils import container_manager
from app.database import models
from app.database.conn import engine, initialize_plugins_from_csv
from celery.signals import worker_shutting_down

def setup_signal_handlers():
    """전역 시그널 핸들러 설정"""
    def signal_handler(signum, frame):
        print(f"Received signal {signum}. Cleaning up all containers...")
        try:
            container_manager.cleanup_all_task_containers()
        except Exception as e:
            print(f"Error during container cleanup: {e}")
        finally:
            sys.exit(0)
    
    # 시그널 핸들러 등록
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

# Celery 이벤트 핸들러를 정의합니다
def on_worker_shut_down(sender=None, conf=None, **kwargs):
    print("Worker is shutting down. Cleaning up containers and trying to reconnect...")
    try:
        # 컨테이너 정리
        container_manager.cleanup_all_task_containers()
        
        # 재연결 시도
        with celery.connection() as connection:
            connection.ensure_connection(max_retries=3)
            print("재연결에 성공했습니다.")
    except Exception as e:
        print(f"Worker shutdown 처리 중 오류 발생: {e}")

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

# 시그널 핸들러 설정
setup_signal_handlers()

def get_celery_app():
    return celery

# 서버 시작 시 CSV 파일로부터 플러그인 데이터를 초기화
@app.on_event("startup")
async def startup_event():
    initialize_plugins_from_csv("./app/database/initial_data/plugin_initialization.csv")