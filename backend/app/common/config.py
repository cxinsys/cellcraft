import secrets
import json
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
        if isinstance(v, str):
            # Check if string starts with '[' - likely JSON array
            if v.startswith("["):
                try:
                    parsed = json.loads(v)
                    if isinstance(parsed, list):
                        return parsed
                    raise ValueError(f"JSON parsed value is not a list: {type(parsed)}")
                except json.JSONDecodeError as e:
                    raise ValueError(f"Invalid JSON format for CORS origins: {e}")
            else:
                # Comma-separated string
                return [i.strip() for i in v.split(",")]
        elif isinstance(v, list):
            return v
        raise ValueError(f"Invalid type for CORS origins: {type(v)}")

    PROJECT_NAME: str = "CellCraft API"
    APP_VERSION: str = "1.0.7"
    APP_DESCRIPTION: str = "Gene Regulatory Network reconstruction platform from scRNA-seq data"
    ENVIRONMENT: str = "local"  # "local" (default) or "production"
    SENTRY_DSN: Optional[HttpUrl] = None
    
    # Pydantic BaseSettings automatically loads these from environment variables
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    POSTGRES_HOST: str
    POSTGRES_PORT: str
    POSTGRES_DB: str

    # Dynamically construct database URI using validator
    SQLALCHEMY_DATABASE_URI: Optional[str] = None

    @validator('SQLALCHEMY_DATABASE_URI', pre=True, always=True)
    def assemble_database_uri(cls, v, values):
        if v:
            return v
        return (
            f"postgresql://{values.get('POSTGRES_USER')}:"
            f"{values.get('POSTGRES_PASSWORD')}@"
            f"{values.get('POSTGRES_HOST')}:"
            f"{values.get('POSTGRES_PORT')}/"
            f"{values.get('POSTGRES_DB')}"
        )

    # Redis Configuration
    REDIS_URL: str = "redis://redis:6379/0"
    DATATABLE_DF_CACHE_TTL: int = 3600  # 1 hour
    PLUGIN_LIST_CACHE_TTL: int = 120    # 2 minutes

    # Resource Queue Configuration
    RESOURCE_CPU_TOTAL: int = 0            # 0 = auto-detect (os.cpu_count()), positive = manual override
    RESOURCE_GPU_TOTAL: int = 0            # 0 = auto-detect (GPU_COUNT env or GPUtil), positive = manual override
    RESOURCE_DEFAULT_CPU_SLOTS: int = 4    # Default CPU slots per task (TENET baseline)
    RESOURCE_DEFAULT_GPU_SLOTS: int = 1    # Default GPU slots per task (FastTENET baseline)
    RESOURCE_POLL_INTERVAL: int = 10       # Retry interval when resources unavailable (seconds)
    CELERY_WORKER_CONCURRENCY: int = 0     # 0 = auto-calculate (cpu_total // default_cpu_slots), positive = manual override

    # Celery configuration with proper defaults
    CELERY_BROKER_URL: str = "amqp://guest:guest@rabbitmq:5672//"
    CELERY_RESULT_BACKEND: Optional[str] = None

    @validator('CELERY_RESULT_BACKEND', pre=True, always=True)
    def assemble_celery_result_backend(cls, v, values):
        if v:
            return v
        db_uri = values.get('SQLALCHEMY_DATABASE_URI')
        if db_uri:
            return f"db+{db_uri}"
        return None

    CELERY_TASK_QUEUES: list = (
        # default queue
        Queue("celery"),
        # custom queue
        Queue("workflow_task"),
        # plugin task queue
        Queue("plugin_task"),
    )

    CELERY_TASK_ROUTES = (route_task,)

    # ==============================================================================
    # File Upload Validation Configuration
    # ==============================================================================

    # Maximum file size: 5GB for H5AD files (large scRNA-seq datasets)
    # Note: sc.read_h5ad() loads entire file into memory during compression,
    # so peak memory ≈ 1.0-1.3x file size. 5GB keeps peak under ~6.5GB.
    MAX_UPLOAD_SIZE: int = 5 * 1024 * 1024 * 1024  # 5GB in bytes

    # Allowed file extensions for upload
    ALLOWED_EXTENSIONS: set = {".h5ad", ".csv", ".json", ".txt"}

    # Maximum filename length (filesystem limit)
    MAX_FILENAME_LENGTH: int = 255

    # Maximum number of files per batch upload request
    MAX_FILES_PER_UPLOAD: int = 20

    # Chunk size for streaming file uploads (1MB chunks)
    # This prevents loading entire files into memory during upload
    # Memory usage: ~1MB per upload instead of up to 500MB
    UPLOAD_CHUNK_SIZE: int = 1024 * 1024  # 1MB in bytes

    # H5AD Compression Configuration
    H5AD_COMPRESSION_ENABLED: bool = True
    H5AD_COMPRESSION_MIN_SIZE: int = 1 * 1024 * 1024  # 1MB

    # ==============================================================================
    # Security Configuration (OWASP Top 10 Compliance)
    # ==============================================================================

    # Security logging configuration (OWASP A09:2021 - Security Logging Failures)
    SECURITY_LOG_FILE: str = "./logs/security.log"

    # Enable file signature validation to prevent file spoofing
    # (OWASP A04:2021 - Insecure Design)
    ENABLE_FILE_SIGNATURE_VALIDATION: bool = True

    # Security logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
    # WARNING: Logs failed attempts and suspicious activity
    # CRITICAL: Logs only critical security events
    SECURITY_LOG_LEVEL: str = "WARNING"

    # Enable anomaly detection in security logs
    # (OWASP A05:2021 - Security Misconfiguration)
    ENABLE_ANOMALY_DETECTION: bool = True

    # Anomaly detection thresholds
    ANOMALY_PATH_TRAVERSAL_THRESHOLD: int = 5  # Max path traversal attempts per hour
    ANOMALY_FILE_SPOOFING_THRESHOLD: int = 3  # Max file signature mismatches per hour
    ANOMALY_USER_EVENTS_THRESHOLD: int = 10  # Max security events per user per hour

    class Config:
        env_file = '.env'
        env_file_encoding = 'utf-8'
        case_sensitive = True

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

# Export file upload validation constants for easy access
MAX_UPLOAD_SIZE = settings.MAX_UPLOAD_SIZE
ALLOWED_EXTENSIONS = settings.ALLOWED_EXTENSIONS
MAX_FILENAME_LENGTH = settings.MAX_FILENAME_LENGTH
MAX_FILES_PER_UPLOAD = settings.MAX_FILES_PER_UPLOAD
UPLOAD_CHUNK_SIZE = settings.UPLOAD_CHUNK_SIZE

# Export security configuration constants
SECURITY_LOG_FILE = settings.SECURITY_LOG_FILE
ENABLE_FILE_SIGNATURE_VALIDATION = settings.ENABLE_FILE_SIGNATURE_VALIDATION
SECURITY_LOG_LEVEL = settings.SECURITY_LOG_LEVEL
ENABLE_ANOMALY_DETECTION = settings.ENABLE_ANOMALY_DETECTION
ANOMALY_PATH_TRAVERSAL_THRESHOLD = settings.ANOMALY_PATH_TRAVERSAL_THRESHOLD
ANOMALY_FILE_SPOOFING_THRESHOLD = settings.ANOMALY_FILE_SPOOFING_THRESHOLD
ANOMALY_USER_EVENTS_THRESHOLD = settings.ANOMALY_USER_EVENTS_THRESHOLD

# Export H5AD compression configuration
H5AD_COMPRESSION_ENABLED = settings.H5AD_COMPRESSION_ENABLED
H5AD_COMPRESSION_MIN_SIZE = settings.H5AD_COMPRESSION_MIN_SIZE