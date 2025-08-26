from typing import List, Dict, Any, Optional, Union
from pydantic import BaseModel
from datetime import datetime
from app.common.enums import PluginType

class DependencyFile(BaseModel):
    file: str  # 파일 내용
    fileName: str  # 파일 이름
    type: str  # 파일 타입

class ReferenceFolder(BaseModel):
    folderName: str  # 폴더 이름
    files: List[DependencyFile]  # 현재 폴더의 파일 리스트
    subFolders: List["ReferenceFolder"] = []  # 하위 폴더 리스트 (재귀적 참조)

    class Config:
        arbitrary_types_allowed = True  # 재귀 참조 허용

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
    pluginType: PluginType  # "analysis" or "visualization"
    referenceFolders: List[ReferenceFolder]
    dependencyFiles: List[DependencyFile] = None
    useGpu: Optional[bool] = False

class PluginData(BaseModel):
    plugin: PluginInfo
    rules: List[Rule]
    drawflow: Dict[str, Any]

class PluginCreate(BaseModel):
    name: str
    description: str
    author: str
    plugin_path: str
    plugin_type: Optional[PluginType] = None  # "analysis" or "visualization"
    dependencies: Optional[Dict[str, str]]
    reference_folders: Optional[Dict[str, Any]] = None
    drawflow: Dict[str, Any]
    rules: Dict[str, Any]
    use_gpu: Optional[bool] = False
    source: Optional[str] = "local"  # "official" or "local"
    is_editable: Optional[bool] = True  # False for official plugins
    version: Optional[str] = None  # Plugin version string
    submodule_path: Optional[str] = None  # Path within submodule for official plugins

class PluginUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    author: Optional[str] = None
    plugin_path: Optional[str] = None
    plugin_type: Optional[PluginType] = None  # "analysis" or "visualization"
    dependencies: Optional[Dict[str, str]] = None
    drawflow: Optional[Dict[str, Any]] = None
    rules: Optional[Dict[str, Any]] = None
    use_gpu: Optional[bool] = None
    source: Optional[str] = None  # "official" or "local"
    is_editable: Optional[bool] = None  # False for official plugins
    version: Optional[str] = None  # Plugin version string
    submodule_path: Optional[str] = None  # Path within submodule for official plugins

class PluginAssociate(BaseModel):
    plugin_id: int

class PluginResponse(BaseModel):
    id: int
    name: str
    description: str
    author: str
    plugin_path: str
    plugin_type: Optional[PluginType] = None  # "analysis" or "visualization"
    dependencies: Optional[Dict[str, str]] = None
    drawflow: Dict[str, Any]
    rules: Dict[str, Any]
    use_gpu: bool = False
    source: str = "local"  # "official" or "local"
    is_editable: bool = True  # False for official plugins
    version: Optional[str] = None  # Plugin version string
    submodule_path: Optional[str] = None  # Path within submodule for official plugins
    created_at: datetime
    updated_at: datetime
    
    class Config:
        from_attributes = True

class BuildDockerRequest(BaseModel):
    use_gpu: Optional[bool] = False
