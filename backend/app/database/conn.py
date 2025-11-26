from contextlib import contextmanager
from datetime import datetime
import pandas as pd
import json
import os
import warnings

from sqlalchemy import Column, DateTime, create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.ext.declarative import declarative_base, as_declarative, declared_attr

from app.common.config import settings
from app.database import models
from app.common.security import get_password_hash

# Use test PostgreSQL database when TESTING=1 (matches production environment)
if os.environ.get("TESTING") == "1":
    # Test DB runs on Docker test-db service (localhost:5433)
    SQLALCHEMY_TEST_DATABASE_URI = "postgresql://test_user:test_pass@localhost:5433/cellcraft_test"
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

@as_declarative()
class Base:
    created_at = Column(DateTime, default=datetime.now)
    updated_at = Column(DateTime, default=datetime.now, onupdate=datetime.now)
    __name__: str

    #CamelCase의 클래스 이름으로부터 snake_case의 테이블 네임 자동 생성
    # @declared_attr
    # def __tablename__(cls) -> str:
    #     return re.sub(r'(?<!^)(?=[A-Z])', '_', cls.__name__).lower()

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

def initialize_plugins_from_csv(csv_file_path: str = None):
    """
    Initialize plugins from CSV file. 
    If no path is provided, reads from official plugins CSV.
    
    Args:
        csv_file_path (str, optional): Path to CSV file. If None, uses official plugins CSV.
        
    Returns:
        int: Number of plugins successfully initialized
    """
    # Default to official plugins CSV if no path provided
    if csv_file_path is None:
        csv_file_path = "./plugin/official/plugins.csv"
    
    # Check if CSV file exists
    if not os.path.exists(csv_file_path):
        print(f"Warning: Plugin CSV file not found at {csv_file_path}")
        return 0
    
    # CSV 파일 읽기
    try:
        df = pd.read_csv(csv_file_path)
    except Exception as e:
        print(f"Error reading CSV file {csv_file_path}: {e}")
        return 0

    # 세션 시작 - get_db_session() context manager 사용
    try:
        with get_db_session() as session:
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
                print("Created admin user")

            # 데이터 추가
            plugins_added = 0
            for index, row in df.iterrows():
                # 플러그인이 이미 존재하는지 확인 (name과 source로 확인)
                source = "official" if "official" in csv_file_path else "local"
                existing_plugin = session.query(models.Plugin).filter_by(
                    name=row['name'],
                    source=source
                ).first()

                # JSON 필드 파싱 (공통 로직)
                dependencies = {}
                if pd.notna(row['dependencies']) and str(row['dependencies']).strip():
                    try:
                        dependencies = json.loads(str(row['dependencies']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in dependencies for plugin {row['name']}, using empty dict")
                        dependencies = {}

                drawflow = {}
                if pd.notna(row['drawflow']) and str(row['drawflow']).strip():
                    try:
                        drawflow = json.loads(str(row['drawflow']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in drawflow for plugin {row['name']}, using empty dict")
                        drawflow = {}

                rules = {}
                if pd.notna(row['rules']) and str(row['rules']).strip():
                    try:
                        rules = json.loads(str(row['rules']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in rules for plugin {row['name']}, using empty dict")
                        rules = {}

                # 추가 필드 파싱 (공통 로직)
                plugin_type = None
                if 'plugin_type' in row and pd.notna(row['plugin_type']) and str(row['plugin_type']).strip():
                    plugin_type = str(row['plugin_type']).strip()

                use_gpu = False
                if 'use_gpu' in row and pd.notna(row['use_gpu']):
                    if isinstance(row['use_gpu'], str):
                        use_gpu = row['use_gpu'].lower().strip() in ('true', '1', 'yes')
                    else:
                        use_gpu = bool(row['use_gpu'])

                if not existing_plugin:
                    try:
                        # Determine plugin path and attributes based on source
                        if source == "official":
                            plugin_path = f"./plugin/official/{row['name']}"
                            is_editable = False
                            # Extract version and submodule path if available
                            version = row.get('version', None) if 'version' in row else None
                            submodule_path = row.get('submodule_path', None) if 'submodule_path' in row else None
                        else:
                            plugin_path = f"./plugin/local/{row['name']}"
                            is_editable = True
                            version = None
                            submodule_path = None

                        # Handle plugin_path from CSV if provided
                        if 'plugin_path' in row and pd.notna(row['plugin_path']) and str(row['plugin_path']).strip():
                            plugin_path = str(row['plugin_path']).strip()

                        # Handle is_editable from CSV if provided (for consistency)
                        if 'is_editable' in row and pd.notna(row['is_editable']):
                            if isinstance(row['is_editable'], str):
                                is_editable = row['is_editable'].lower().strip() in ('true', '1', 'yes')
                            else:
                                is_editable = bool(row['is_editable'])

                        # Handle author field properly
                        author = str(row['author']) if pd.notna(row['author']) else 'Unknown'

                        plugin = models.Plugin(
                            name=str(row['name']),
                            description=str(row['description']) if pd.notna(row['description']) else '',
                            author=author,
                            plugin_path=plugin_path,
                            plugin_type=plugin_type,
                            dependencies=dependencies,
                            drawflow=drawflow,
                            rules=rules,
                            use_gpu=use_gpu,
                            source=source,
                            is_editable=is_editable,
                            version=version,
                            submodule_path=submodule_path
                        )
                        session.add(plugin)
                        plugins_added += 1
                        print(f"Added {source} plugin: {row['name']} (type: {plugin_type}, GPU: {use_gpu})")
                    except Exception as e:
                        print(f"Error creating plugin {row['name']}: {e}")
                        continue
                else:
                    # Update existing plugin with CSV data
                    try:
                        updated = False

                        # Update fields that might have changed
                        if existing_plugin.description != str(row['description']):
                            existing_plugin.description = str(row['description']) if pd.notna(row['description']) else ''
                            updated = True

                        if existing_plugin.plugin_type != plugin_type:
                            existing_plugin.plugin_type = plugin_type
                            updated = True

                        if existing_plugin.use_gpu != use_gpu:
                            existing_plugin.use_gpu = use_gpu
                            updated = True

                        # Update JSON fields if they differ
                        if existing_plugin.dependencies != dependencies:
                            existing_plugin.dependencies = dependencies
                            updated = True

                        if existing_plugin.drawflow != drawflow:
                            existing_plugin.drawflow = drawflow
                            updated = True

                        if existing_plugin.rules != rules:
                            existing_plugin.rules = rules
                            updated = True

                        # Update version info for official plugins
                        if source == "official":
                            version = row.get('version', None) if 'version' in row else None
                            submodule_path = row.get('submodule_path', None) if 'submodule_path' in row else None

                            if existing_plugin.version != version:
                                existing_plugin.version = version
                                updated = True

                            if existing_plugin.submodule_path != submodule_path:
                                existing_plugin.submodule_path = submodule_path
                                updated = True

                        if updated:
                            print(f"Updated {source} plugin: {row['name']} (type: {plugin_type}, GPU: {use_gpu})")
                        else:
                            print(f"Plugin {row['name']} from {source} is up to date")

                    except Exception as e:
                        print(f"Error updating plugin {row['name']}: {e}")
                        continue

            # 커밋은 context manager가 자동 처리
            print(f"Successfully initialized {plugins_added} plugins from {csv_file_path}")
            return plugins_added
    except Exception as e:
        # rollback/close는 context manager가 자동 처리
        print(f"Error during plugin initialization: {e}")
        return 0