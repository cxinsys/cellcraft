from typing import Optional, Dict, Any, List
from pydantic import BaseModel
from datetime import datetime


class PluginInfo(BaseModel):
    """Plugin information for task monitoring response."""
    version: Optional[str] = None
    source: Optional[str] = None
    plugin_type: Optional[str] = None


class TaskMonitoringItem(BaseModel):
    """Individual task item in task monitoring response."""
    id: int
    task_id: str
    workflow_id: Optional[int] = None  # Allow None for tasks without associated workflow
    user_id: int
    status: str
    start_time: datetime
    end_time: Optional[datetime] = None
    workflow_title: Optional[str] = None
    task_title: Optional[str] = None
    plugin_name: Optional[str] = None
    task_type: Optional[str] = None
    plugin: Optional[PluginInfo] = None


class TaskMonitoringResponse(BaseModel):
    """Response schema for task monitoring API endpoint."""
    tasks: List[TaskMonitoringItem]
    
    class Config:
        orm_mode = True