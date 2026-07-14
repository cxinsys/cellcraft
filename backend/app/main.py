from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware

from app.api import api_router
from app.core.config import settings
from app.core import startup
from app.worker.app import celery


def get_celery_app():
    return celery


def create_app() -> FastAPI:
    startup.run_migrations()

    is_production = settings.ENVIRONMENT == "production"

    app = FastAPI(
        title=settings.PROJECT_NAME,
        version=settings.APP_VERSION,
        description=settings.APP_DESCRIPTION,
        docs_url=None if is_production else "/docs",
        redoc_url=None if is_production else "/redoc",
        openapi_url=None if is_production else "/openapi.json",
    )

    cors_origins = (
        ["https://cellcraft.app"]
        if is_production
        else [str(origin) for origin in settings.BACKEND_CORS_ORIGINS]
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=cors_origins,
        allow_credentials=True,
        allow_methods=["GET", "POST", "PUT", "DELETE", "PATCH"],
        allow_headers=["Authorization", "Content-Type"],
    )

    app.celery_app = celery
    app.include_router(api_router, prefix=settings.ROUTES_STR)
    app.add_event_handler("startup", startup.on_startup)
    startup.setup_signal_handlers()

    return app


app = create_app()
