from sqlalchemy import Table, Boolean, Column, ForeignKey, Integer, String, Text, DateTime, UniqueConstraint, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB, JSON, ARRAY
from sqlalchemy.sql import func

from app.database.conn import Base
from app.common.enums import PluginType

user_plugin_association = Table(
    'user_plugin_association', Base.metadata,
    Column('user_id', Integer, ForeignKey('users.id')),
    Column('plugin_id', Integer, ForeignKey('plugins.id'))
)

class User(Base):
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    username = Column(String, nullable=False)
    email = Column(String, unique=True, index=True, nullable=False)
    hashed_password = Column(String, nullable=False)
    is_active = Column(Boolean, default=True)
    is_superuser = Column(Boolean, default=False)

    files = relationship("File", back_populates="user")
    workflows = relationship("Workflow", back_populates="user")
    tasks = relationship("Task", back_populates="user")
    plugins = relationship("Plugin", secondary=user_plugin_association, back_populates="users")

class File(Base):
    __tablename__ = "files"

    id = Column(Integer, primary_key=True, index=True)
    file_name = Column(String, nullable=False)
    file_size = Column(String, nullable=False)
    file_path = Column(String, nullable=False)
    folder = Column(String, nullable=False)
    user_id = Column(Integer, ForeignKey("users.id"))

    user = relationship("User", back_populates="files")

class Workflow(Base):
    __tablename__ = "workflows"
    # __table_args__ = (UniqueConstraint('user_id', 'title', name='uix_1'), )


    id = Column(Integer, primary_key=True, index=True)
    title = Column(String, nullable=False)
    thumbnail = Column(Text, nullable=True)
    workflow_info = Column(JSONB, nullable=False)
    user_id = Column(Integer, ForeignKey("users.id"))

    user = relationship("User", back_populates="workflows")
    tasks = relationship("Task", back_populates="workflows")

class Task(Base):
    __tablename__ = 'tasks'

    id = Column(Integer, primary_key=True)
    task_id = Column(String, nullable=False)
    start_time = Column(DateTime(timezone=True), nullable=False)
    end_time = Column(DateTime(timezone=True), nullable=True)
    status = Column(String, nullable=False)
    user_id = Column(Integer, ForeignKey("users.id"))
    workflow_id = Column(Integer, ForeignKey("workflows.id"))
    algorithm_id = Column(String, nullable=True)  # algorithm_id 필드 추가
    plugin_name = Column(String, nullable=True)  # plugin_name 필드 추가 (kept for backward compatibility)
    plugin_id = Column(Integer, ForeignKey("plugins.id"), nullable=True)  # New FK to plugins table
    task_type = Column(String, nullable=True)  # task_type 필드 추가 (compile, visualization)
    plugin_image_uri = Column(String, nullable=True)  # Plugin image URI for reproducibility

    user = relationship("User", back_populates="tasks")
    workflows = relationship("Workflow", back_populates="tasks")
    plugin = relationship("Plugin", back_populates="tasks")

class Plugin(Base):
    __tablename__ = "plugins"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, nullable=False)
    description = Column(String, nullable=False)
    author = Column(String, nullable=False)
    plugin_path = Column(String, nullable=False)
    plugin_type = Column(Enum(PluginType), nullable=True)  # "analysis" or "visualization"
    dependencies = Column(JSONB, nullable=True)
    drawflow = Column(JSONB, nullable=False)
    rules = Column(JSONB, nullable=False)
    use_gpu = Column(Boolean, default=False)
    source = Column(String, default="local", nullable=False)  # "official" or "local"
    is_editable = Column(Boolean, default=True, nullable=False)  # False for official plugins
    version = Column(String, nullable=True)  # Plugin version string
    submodule_path = Column(String, nullable=True)  # Path within submodule for official plugins
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    users = relationship("User", secondary=user_plugin_association, back_populates="plugins")
    tasks = relationship("Task", back_populates="plugin")
