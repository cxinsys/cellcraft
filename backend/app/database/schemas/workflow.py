from typing import Dict, List, Optional
from pydantic import BaseModel

class WorkflowBase(BaseModel):
    title: Optional[str] = None
    thumbnail: Optional[str] = None
    workflow_info: Optional[Dict] = None 

class WorkflowCreate(WorkflowBase):
    id: int = None
    title: str = None
    thumbnail: str = None
    workflow_info: Dict = None
    option_file: str = None

class WorkflowFind(WorkflowBase):
    id: int

class WorkflowDelete(WorkflowBase):
    id: int

class Workflow(WorkflowBase):
    id: int
    user_id: int
    
    class Config:
        orm_mode = True

class WorkflowResult(BaseModel):
    filename: Optional[str] = None

class WorkflowNodeFileBase(BaseModel):
    id: Optional[int] = None
    node_id: Optional[str] = None
    node_name: Optional[str] = None
    file_content: Optional[Dict] = None
    file_extension: Optional[str] = None

class WorkflowNodeFileCreate(WorkflowNodeFileBase):
    id: int
    node_id: str
    node_name: str
    file_content: Optional[Dict]
    file_extension: str

class WorkflowNodeFileDelete(WorkflowNodeFileBase):
    id: int
    node_id: str

class WorkflowNodeFileUpdate(WorkflowNodeFileBase):
    id: int
    node_id: str
    file_content: Optional[Dict]

class WorkflowNodeFileRead(WorkflowNodeFileBase):
    id: int
    node_id: str