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

class FileGet(FileBase):
    file_name: str
    local_file_path: str

class FileFind(FileBase):
    file_name: str
    anno_column: Optional[str] = None

class FolderFind(FileBase):
    folder_name: str

class FileDelete(FileBase):
    file_name: str

class FileUpdate(FileBase):
    file_name: str
    update_name: str

class FileResultFind(FileBase):
    file_name: str
    option_file_name: str

class File(FileBase):
    id: int
    user_id: int
    
    class Config:
        orm_mode = True
