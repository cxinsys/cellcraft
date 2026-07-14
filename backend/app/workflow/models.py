from sqlalchemy import Column, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB

from app.db.base import Base


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
