from datetime import datetime
import pandas as pd
from sqlalchemy import Column, DateTime, create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.ext.declarative import declarative_base, as_declarative, declared_attr

from app.common.config import settings
from app.database import models
from app.common.security import get_password_hash

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

def initialize_plugins_from_csv(csv_file_path: str):
    # CSV 파일 읽기
    df = pd.read_csv(csv_file_path)
    
    # 세션 시작
    session = get_new_engine_and_session()

    # 데이터 추가
    for index, row in df.iterrows():
        # 플러그인이 이미 존재하는지 확인
        existing_plugin = session.query(models.Plugin).filter_by(name=row['name']).first()
        if not existing_plugin:
            plugin = models.Plugin(
                name=row['name'],
                description=row['description'],
                author=row['author'],
                plugin_path=row['plugin_path'],
                dependencies=row['dependencies'],  # 필요에 따라 변환 필요
                drawflow=row['drawflow'],  # 필요에 따라 변환 필요
                rules=row['rules']  # 필요에 따라 변환 필요
            )
            session.add(plugin)

    # 관리자 사용자 추가
    existing_user = session.query(models.User).filter_by(username="admin").first()
    if not existing_user:
        hashed_password = get_password_hash("cellcraft2024!")
        user = models.User(
            username="admin",
            email="cellcraft@cellcraft.com",
            hashed_password=hashed_password,
            is_active=True,
            is_superuser=True
        )
        session.add(user)

    # 커밋하여 변경 사항 저장
    session.commit()

    # 세션 종료
    session.close()