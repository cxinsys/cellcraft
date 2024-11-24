from typing import List, Dict, Any, Optional, Union
from pydantic import BaseModel

class DependencyFile(BaseModel):
    file: str
    fileName: str
    type: str

class Parameter(BaseModel):
    name: str
    type: str
    defaultValue: Optional[Any] = None
    min: Optional[Any] = None
    max: Optional[Any] = None
    fileExtension: Optional[str] = None

class Rule(BaseModel):
    name: str
    input: List[str]
    output: List[str]
    script: Optional[str] = None
    parameters: List[Parameter]
    nodeId: int
    isVisualization: Optional[bool] = False

class PluginInfo(BaseModel):
    name: str
    description: str
    dependencyFiles: List[DependencyFile]

class PluginData(BaseModel):
    plugin: PluginInfo
    rules: List[Rule]
    drawflow: Dict[str, Any]

class PluginCreate(BaseModel):
    name: str
    description: str
    author: str
    plugin_path: str
    dependencies: Optional[Dict[str, str]]
    drawflow: Dict[str, Any]
    rules: Dict[str, Any]

class PluginUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    author: Optional[str] = None
    plugin_path: Optional[str] = None
    dependencies: Optional[Dict[str, str]] = None
    drawflow: Optional[Dict[str, Any]] = None
    rules: Optional[Dict[str, Any]] = None

class PluginAssociate(BaseModel):
    plugin_id: int
