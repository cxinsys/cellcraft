from typing import Dict, List, Optional
from typing import Optional
from pydantic import BaseModel, EmailStr

class WorkflowBase(BaseModel):
    title: Optional[str] = None
    workflow_info: Optional[Dict[str, float]] = None 

class WorkflowCreate(WorkflowBase):
    title: str
    workflow_info: Dict[str, float] = None

class Workflow(WorkflowBase):
    id: int
    user_id: int
    
    class Config:
        orm_mode = True

class WorkflowResult(BaseModel):
    filename: Optional[str] = None