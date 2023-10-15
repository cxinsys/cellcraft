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

class FileSetup(FileBase):
    file_name: str
    anno_of_interest: str
    pseudo_of_interest: str
    clusters_of_interest: List[str]
    make_binary: Optional[bool] = None
    device: Optional[str] = None
    device_ids: Optional[List] = None
    batch_size: Optional[int] = None
    kp: Optional[float] = None
    percentile: Optional[int] = None
    win_length: Optional[int] = None
    polyorder: Optional[int] = None

class File(FileBase):
    id: int
    user_id: int
    
    class Config:
        orm_mode = True
