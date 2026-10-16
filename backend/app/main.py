from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from starlette.middleware.cors import CORSMiddleware

from app.api import api_router
from app.core.config import settings
from app.core import startup
from app.core.exceptions import CellcraftError
from app.worker.app import celery


def get_celery_app():
    return celery


async def cellcraft_error_handler(request: Request, exc: CellcraftError) -> JSONResponse:
    """Map a domain ``CellcraftError`` onto FastAPI's HTTPException wire format.

    Emits ``{"detail": <detail>}`` with the exception's status code — byte-for-byte
    identical to FastAPI's default ``HTTPException`` response — so converting a
    service ``raise HTTPException(...)`` to the matching domain exception does not
    change any response. FastAPI's built-in 422 ``RequestValidationError`` handling
    is untouched (this handler only fires for ``CellcraftError`` subclasses).
    """
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.detail})


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
    app.add_exception_handler(CellcraftError, cellcraft_error_handler)
    app.include_router(api_router, prefix=settings.ROUTES_STR)
    app.add_event_handler("startup", startup.on_startup)
    startup.setup_signal_handlers()

    return app


app = create_app()
