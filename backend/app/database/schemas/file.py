from typing import List, Optional
from typing import Optional
from pydantic import BaseModel, EmailStr

class FileBase(BaseModel):
    file_name: Optional[str] = None
    file_size: Optional[str] = None
    file_path: Optional[str] = None

class FileCreate(FileBase):
    file_name: str
    file_size: str
    file_path: str

class File(FileBase):
    id: int
    user_id: int
    
    class Config:
        orm_mode = True
