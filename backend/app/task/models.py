from sqlalchemy import Column, ForeignKey, Integer, String, DateTime
from sqlalchemy.orm import relationship

from app.db.base import Base


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
    # index=True matches ix_tasks_plugin_id created by migration 0001
    # (historically missing from the model).
    plugin_id = Column(Integer, ForeignKey("plugins.id"), nullable=True, index=True)
    task_type = Column(String, nullable=True)  # task_type 필드 추가 (compile, visualization)
    plugin_image_uri = Column(String, nullable=True)  # Plugin image URI for reproducibility

    user = relationship("User", back_populates="tasks")
    workflows = relationship("Workflow", back_populates="tasks")
    plugin = relationship("Plugin", back_populates="tasks")
