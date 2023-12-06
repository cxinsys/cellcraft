from typing import Dict, List, Optional
from typing import Optional
from pydantic import BaseModel, EmailStr

class WorkflowBase(BaseModel):
    title: Optional[str] = None
    thumbnail: Optional[str] = None
    workflow_info: Optional[Dict] = None 

class WorkflowCreate(WorkflowBase):
    id: int = None
    title: str = None
    thumbnail: str = None
    workflow_info: Dict = None
    nodes: List = None
    linked_nodes: List = None

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