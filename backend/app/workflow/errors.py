"""Workflow 유틸리티 예외 계층 (PR-10 split from ``workflow/utils``).

기존 ``workflow/utils.py`` 상단에 정의돼 있던 예외들을 별도 모듈로 분리한다.
Snakefile 처리 함수들이 ``compiler/visualization.py``로 이동하면서, 이동한 함수와
잔류한 함수(``load_tab_file`` 등)가 동일한 예외를 공유해야 하는데, utils 파사드가
compiler를 re-export 하므로 예외를 utils에 두면 순환 import가 생긴다. 이를 피하기
위해 예외를 leaf 모듈로 분리했다. ``app.workflow.utils``가 그대로 re-export 하므로
기존 참조 경로(``app.workflow.utils.FileValidationError`` 등)는 유지된다.
"""
import logging
from typing import Dict, List, Any

# Configure logging
logger = logging.getLogger(__name__)


# Custom exceptions for better error handling
class WorkflowUtilsError(Exception):
    """Base exception for workflow utilities."""
    pass

class FileValidationError(WorkflowUtilsError):
    """Exception raised when file validation fails."""
    def __init__(self, message: str, file_path: str = None, context: Dict[str, Any] = None):
        super().__init__(message)
        self.file_path = file_path
        self.context = context or {}
        logger.error(f"File validation error: {message}. Path: {file_path}. Context: {context}")

class SnakefileProcessingError(WorkflowUtilsError):
    """Exception raised when Snakefile processing fails."""
    def __init__(self, message: str, snakefile_path: str = None, rule_name: str = None):
        super().__init__(message)
        self.snakefile_path = snakefile_path
        self.rule_name = rule_name
        logger.error(f"Snakefile processing error: {message}. File: {snakefile_path}, Rule: {rule_name}")

class AlgorithmResolutionError(WorkflowUtilsError):
    """Exception raised when algorithm path resolution fails."""
    def __init__(self, message: str, node_id: str = None, available_files: List = None):
        super().__init__(message)
        self.node_id = node_id
        self.available_files = available_files
        logger.error(f"Algorithm resolution error: {message}. Node ID: {node_id}")
