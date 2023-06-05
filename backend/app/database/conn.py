from datetime import datetime
import re
from click import echo
from sqlalchemy import Column, DateTime, create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.ext.declarative import declarative_base, as_declarative, declared_attr

from app.common.config import settings

engine = create_engine(settings.SQLALCHEMY_DATABASE_URI, echo=True ,pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

@as_declarative()
class Base:
    created_at = Column(DateTime, default=datetime.now)
    updated_at = Column(DateTime, default=datetime.now, onupdate=datetime.now)
    __name__: str

    #CamelCase의 클래스 이름으로부터 snake_case의 테이블 네임 자동 생성
    # @declared_attr
    # def __tablename__(cls) -> str:
    #     return re.sub(r'(?<!^)(?=[A-Z])', '_', cls.__name__).lower()

def get_new_engine_and_session() -> Session:
    engine = create_engine(settings.SQLALCHEMY_DATABASE_URI, pool_pre_ping=True)
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    return SessionLocal()