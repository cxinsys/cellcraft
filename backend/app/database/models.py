from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, func
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSON


from app.database.conn import Base


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

class File(Base):
    __tablename__ = "files"

    id = Column(Integer, primary_key=True, index=True)
    file_name = Column(String, nullable=False)
    file_size = Column(String, nullable=False)
    file_path = Column(String, nullable=False)
    user_id = Column(Integer, ForeignKey("users.id"))

    user = relationship("User", back_populates="files")

class Workflow(Base):
    __tablename__ = "workflows"

    id = Column(Integer, primary_key=True, index=True)
    title = Column(String, nullable=False)
    workflow_info = Column(JSON, nullable=False)
    user_id = Column(Integer, ForeignKey("users.id"))

    user = relationship("User", back_populates="workflows")