from typing import List, Dict, Any, Optional, Union
from pydantic import BaseModel

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
    referenceFolders: List[ReferenceFolder]
    dependencyFiles: List[DependencyFile] = None

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
    reference_folders: Optional[Dict[str, Union[str, Dict[str, Any]]]] = None
    drawflow: Dict[str, Any]
    rules: Dict[str, Any]

class PluginUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    author: Optional[str] = None
    plugin_path: Optional[str] = None
    dependencies: Optional[Dict[str, str]] = None
    reference_folders: Optional[Dict[str, Union[str, Dict[str, Any]]]] = None
    drawflow: Optional[Dict[str, Any]] = None
    rules: Optional[Dict[str, Any]] = None

class PluginAssociate(BaseModel):
    plugin_id: int
