import secrets
from typing import Any, Dict, List, Optional, Union
from os import environ
from dotenv import load_dotenv
from functools import lru_cache
from kombu import Queue
from pydantic import AnyHttpUrl, BaseSettings, EmailStr, HttpUrl, PostgresDsn, validator

load_dotenv(verbose=True)

def route_task(name, args, kwargs, options, task=None, **kw):
    if ":" in name:
        queue, _ = name.split(":")
        return {"queue": queue}
    return {"queue": "celery"}

class Settings(BaseSettings):
    ROUTES_STR: str = "/routes"
    SECRET_KEY: str = secrets.token_urlsafe(32)
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24 * 8
    SERVER_NAME: str = None
    SERVER_HOST: AnyHttpUrl = None

    BACKEND_CORS_ORIGINS: List[AnyHttpUrl] = []

    @validator("BACKEND_CORS_ORIGINS", pre=True)
    def assemble_cors_origins(cls, v: Union[str, List[str]]) -> Union[List[str], str]:
        if isinstance(v, str) and not v.startswith("["):
            return [i.strip() for i in v.split(",")]
        elif isinstance(v, (list, str)):
            return v
        raise ValueError(v)

    PROJECT_NAME: str = "test"
    SENTRY_DSN: Optional[HttpUrl] = None
    
    POSTGRES_USER: str = environ.get('POSTGRES_USER')
    POSTGRES_PASSWORD: str = environ.get('POSTGRES_PASSWORD')
    POSTGRES_HOST: str = environ.get('POSTGRES_HOST')
    POSTGRES_PORT: str = environ.get('POSTGRES_PORT')
    POSTGRES_DB: str = environ.get('POSTGRES_DB')

    SQLALCHEMY_DATABASE_URI: str = f"postgresql://{POSTGRES_USER}:{POSTGRES_PASSWORD}@{POSTGRES_HOST}:{POSTGRES_PORT}/{POSTGRES_DB}"

    CELERY_BROKER_URL: str = environ.get("CELERY_BROKER_URL", "amqp://guest:guest@rabbitmq:5672//")
    CELERY_RESULT_BACKEND: str = environ.get("CELERY_RESULT_BACKEND", f"db+{SQLALCHEMY_DATABASE_URI}")

    CELERY_TASK_QUEUES: list = (
        # default queue
        Queue("celery"),
        # custom queue
        Queue("workflow_task"),
    )

    CELERY_TASK_ROUTES = (route_task,)

class DevelopmentConfig(Settings):
    pass

@lru_cache()
def get_settings():
    config_cls_dict = {
        "development": DevelopmentConfig,
    }
    config_name = environ.get("CELERY_CONFIG", "development")
    config_cls = config_cls_dict[config_name]
    return config_cls()

settings = get_settings()