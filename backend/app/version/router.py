from fastapi import APIRouter

from app.core.config import settings

router = APIRouter()


@router.get("")
def get_version():
    """애플리케이션 버전 정보를 반환한다. 인증 불필요."""
    return {
        "version": settings.APP_VERSION,
        "name": settings.PROJECT_NAME,
        "environment": settings.ENVIRONMENT,
    }
