from sqlalchemy import Table, Boolean, Column, ForeignKey, Integer, String, false
from sqlalchemy.orm import relationship

from app.db.base import Base

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
    # 초기·임시 비밀번호 사용 중 → 비밀번호 변경 전까지 다른 API 사용 불가 (#138)
    must_change_password = Column(Boolean, nullable=False, default=False, server_default=false())

    files = relationship("File", back_populates="user")
    workflows = relationship("Workflow", back_populates="user")
    tasks = relationship("Task", back_populates="user")
    plugins = relationship("Plugin", secondary=user_plugin_association, back_populates="users")
