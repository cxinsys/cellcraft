from sqlalchemy import Boolean, Column, String, Integer, DateTime, Enum, UniqueConstraint
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.sql import func

from app.db.base import Base
from app.core.enums import PluginType
from app.user.models import user_plugin_association


class Plugin(Base):
    __tablename__ = "plugins"
    # Declared in migration 0001 but historically missing from the model —
    # declared here so autogenerate sees model == deployed schema.
    __table_args__ = (
        UniqueConstraint('name', 'source', name='uq_plugin_name_source'),
    )

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
