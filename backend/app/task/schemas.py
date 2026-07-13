from typing import Optional, Dict, Any, List
from pydantic import BaseModel
from datetime import datetime


class TaskSubmission(BaseModel):
    """Typed representation of the ``process_data_task`` Celery dispatch kwargs.

    Extracted in PR-7 (Phase 3c). This model types the keyword payload sent to
    ``app.worker.tasks.process_data_task`` so the dispatch site can build the
    ``kwargs`` dict from a validated schema instead of an inline literal.

    Wire-format contract (IMPORTANT): the on-the-wire ``kwargs`` keys, their
    names, and their value shapes are preserved exactly for old-worker /
    new-web compatibility. ``resource_type`` / ``resource_slots`` carry the same
    defaults as the worker signature (``'cpu'`` / ``4``). ``cache_key`` and
    ``cache_info`` are the visualization-only extras absorbed by the worker's
    ``**kwargs``; they are ``None`` for compile tasks. Callers therefore emit
    the exact legacy dict via ``submission.dict(exclude_none=True)``:

    - compile:       user_id, workflow_id, algorithm_id, plugin_name,
                     task_type, resource_type, resource_slots
    - visualization: the above + cache_key, cache_info

    pydantic v1 syntax only.
    """
    user_id: int
    workflow_id: int
    algorithm_id: Any
    plugin_name: str
    task_type: str
    resource_type: str = 'cpu'
    resource_slots: int = 4
    cache_key: Optional[str] = None
    cache_info: Optional[Dict[str, Any]] = None


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