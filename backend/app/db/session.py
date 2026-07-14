from contextlib import contextmanager
from typing import Generator
import os
import warnings

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session

from app.core.config import settings
from app.db.base import Base  # noqa: F401 — re-exported for backward-compatible imports

# Use test PostgreSQL database when TESTING=1 (matches production environment)
if os.environ.get("TESTING") == "1":
    # Test DB runs on Docker test-db service (localhost:5433 by default).
    # TEST_DATABASE_URI overrides it when tests run inside a container
    # (e.g. docker-compose.test.yml uses test-db:5432 on its own network).
    SQLALCHEMY_TEST_DATABASE_URI = os.environ.get(
        "TEST_DATABASE_URI",
        "postgresql://test_user:test_pass@localhost:5433/cellcraft_test",
    )
    engine = create_engine(
        SQLALCHEMY_TEST_DATABASE_URI,
        echo=False,  # Disable SQL logging in tests
        pool_pre_ping=True,
        pool_recycle=3600  # Recycle connections after 1 hour
    )
else:
    engine = create_engine(
        settings.SQLALCHEMY_DATABASE_URI,
        echo=True,
        pool_pre_ping=True,
        pool_recycle=3600  # Recycle connections after 1 hour
    )

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


# FastAPI request-scoped DB dependency (relocated from the legacy routes.dep in PR-3).
def get_db() -> Generator:
    try:
        db = SessionLocal()
        yield db
    finally:
        db.close()


@contextmanager
def get_db_session():
    """
    Database session context manager - 글로벌 Engine 재사용.
    Celery 태스크 및 백그라운드 작업용.

    Usage:
        with get_db_session() as db:
            db.query(Model).all()
            # commit은 자동으로 처리됨
    """
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def get_new_engine_and_session() -> Session:
    """
    DEPRECATED: get_db_session() context manager를 사용하세요.
    이 함수는 Engine 누수를 발생시킵니다.

    하위 호환성을 위해 유지되지만, 글로벌 Engine을 재사용하도록 변경되었습니다.
    """
    warnings.warn(
        "get_new_engine_and_session()은 deprecated입니다. "
        "get_db_session() context manager를 사용하세요.",
        DeprecationWarning,
        stacklevel=2
    )
    return SessionLocal()  # 글로벌 엔진 사용으로 변경

