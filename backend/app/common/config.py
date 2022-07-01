import secrets
from typing import Any, Dict, List, Optional, Union
from os import environ
from dotenv import load_dotenv

from pydantic import AnyHttpUrl, BaseSettings, EmailStr, HttpUrl, PostgresDsn, validator

load_dotenv(verbose=True)

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
    
    POSTGRES_USERNAME: str = environ.get('POSTGRES_USERNAME')
    POSTGRES_PASSWORD: str = environ.get('POSTGRES_PASSWORD')
    POSTGRES_HOST: str = environ.get('POSTGRES_HOST')
    POSTGRES_PORT: str = environ.get('POSTGRES_PORT')
    POSTGRES_DB: str = environ.get('POSTGRES_DB')

    SQLALCHEMY_DATABASE_URI: str = f"postgresql://{POSTGRES_USERNAME}:{POSTGRES_PASSWORD}@{POSTGRES_HOST}:{POSTGRES_PORT}/{POSTGRES_DB}"

settings = Settings()