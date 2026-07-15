"""
Error handling utilities for CellCraft backend API.

This module provides consistent error response formats and specialized HTTP exceptions
for better error handling across the application.
"""

from typing import Any, Dict, Optional, List
from fastapi import HTTPException, status
from fastapi.responses import JSONResponse
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ErrorCategory(str, Enum):
    """Error categories for consistent error classification."""
    VALIDATION = "validation"
    NOT_FOUND = "not_found"
    PERMISSION = "permission"
    SERVER_ERROR = "server_error"
    PLUGIN_ERROR = "plugin_error"
    FILE_ERROR = "file_error"
    WORKFLOW_ERROR = "workflow_error"
    TASK_ERROR = "task_error"


class ErrorResponse:
    """Consistent error response format for API endpoints."""
    
    def __init__(
        self,
        error_type: ErrorCategory,
        message: str,
        details: Optional[str] = None,
        suggested_actions: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None
    ):
        self.error_type = error_type
        self.message = message
        self.details = details
        self.suggested_actions = suggested_actions or []
        self.context = context or {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert error response to dictionary format."""
        response = {
            "error": {
                "type": self.error_type.value,
                "message": self.message,
                "details": self.details,
                "suggested_actions": self.suggested_actions,
                "context": self.context
            }
        }
        return response


class CellCraftHTTPException(HTTPException):
    """Enhanced HTTP exception with structured error responses."""
    
    def __init__(
        self,
        status_code: int,
        error_type: ErrorCategory,
        message: str,
        details: Optional[str] = None,
        suggested_actions: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None
    ):
        self.error_response = ErrorResponse(
            error_type=error_type,
            message=message,
            details=details,
            suggested_actions=suggested_actions,
            context=context
        )
        super().__init__(status_code=status_code, detail=self.error_response.to_dict())


# Specialized Exception Classes

class PluginNotFoundError(CellCraftHTTPException):
    """Plugin not found error."""
    
    def __init__(self, plugin_name: str, available_plugins: Optional[List[str]] = None):
        suggested_actions = [
            "Check if the plugin name is spelled correctly",
            "Verify the plugin is installed and available",
            "Contact administrator if plugin should be available"
        ]
        
        if available_plugins:
            suggested_actions.append(f"Available plugins: {', '.join(available_plugins[:5])}")
        
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            error_type=ErrorCategory.NOT_FOUND,
            message=f"Plugin '{plugin_name}' not found",
            details=f"The requested plugin '{plugin_name}' is not available in the system",
            suggested_actions=suggested_actions,
            context={"plugin_name": plugin_name, "available_plugins": available_plugins}
        )


class ScriptNotFoundError(CellCraftHTTPException):
    """Visualization script not found error."""
    
    def __init__(self, script_name: str, plugin_name: str, available_scripts: Optional[List[str]] = None):
        suggested_actions = [
            "Check if the script name is correct",
            f"Verify the script exists in plugin '{plugin_name}'",
            "Try selecting a different visualization script"
        ]
        
        if available_scripts:
            suggested_actions.append(f"Available scripts: {', '.join(available_scripts[:5])}")
        
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            error_type=ErrorCategory.NOT_FOUND,
            message=f"Script '{script_name}' not found in plugin '{plugin_name}'",
            details=f"The visualization script '{script_name}' does not exist in the specified plugin",
            suggested_actions=suggested_actions,
            context={
                "script_name": script_name,
                "plugin_name": plugin_name,
                "available_scripts": available_scripts
            }
        )


class FileNotFoundError(CellCraftHTTPException):
    """File not found error."""
    
    def __init__(self, file_path: str, file_type: str = "file"):
        suggested_actions = [
            "Check if the file path is correct",
            "Verify the file exists and is accessible",
            "Ensure the workflow has been executed successfully",
            "Check file permissions"
        ]
        
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            error_type=ErrorCategory.FILE_ERROR,
            message=f"{file_type.title()} not found",
            details=f"The requested {file_type} '{file_path}' does not exist or is not accessible",
            suggested_actions=suggested_actions,
            context={"file_path": file_path, "file_type": file_type}
        )


class ValidationError(CellCraftHTTPException):
    """Input validation error."""
    
    def __init__(self, field_name: str, message: str, required_format: Optional[str] = None):
        suggested_actions = [
            f"Check the '{field_name}' field",
            "Ensure all required fields are provided",
            "Verify the data format is correct"
        ]
        
        if required_format:
            suggested_actions.append(f"Required format: {required_format}")
        
        super().__init__(
            status_code=status.HTTP_400_BAD_REQUEST,
            error_type=ErrorCategory.VALIDATION,
            message=f"Invalid {field_name}",
            details=message,
            suggested_actions=suggested_actions,
            context={"field_name": field_name, "required_format": required_format}
        )


class WorkflowError(CellCraftHTTPException):
    """Workflow-related error."""
    
    def __init__(self, workflow_id: str, message: str, error_type: str = "workflow"):
        suggested_actions = [
            "Check workflow configuration",
            "Verify workflow exists and is accessible",
            "Ensure proper workflow permissions"
        ]
        
        super().__init__(
            status_code=status.HTTP_400_BAD_REQUEST,
            error_type=ErrorCategory.WORKFLOW_ERROR,
            message=f"Workflow error: {error_type}",
            details=message,
            suggested_actions=suggested_actions,
            context={"workflow_id": workflow_id, "error_type": error_type}
        )


class SnakefileGenerationError(CellCraftHTTPException):
    """Snakefile generation error."""
    
    def __init__(self, plugin_name: str, script_name: str, error_details: str):
        suggested_actions = [
            "Check plugin configuration",
            "Verify script parameters are valid",
            "Contact administrator if problem persists",
            "Try with different parameters"
        ]
        
        super().__init__(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_type=ErrorCategory.SERVER_ERROR,
            message="Failed to generate visualization Snakefile",
            details=f"Could not generate Snakefile for '{script_name}' in plugin '{plugin_name}': {error_details}",
            suggested_actions=suggested_actions,
            context={
                "plugin_name": plugin_name,
                "script_name": script_name,
                "error_details": error_details
            }
        )


class TaskSubmissionError(CellCraftHTTPException):
    """Task submission error."""
    
    def __init__(self, task_type: str, error_details: str):
        suggested_actions = [
            "Check system resources",
            "Try again in a few minutes",
            "Contact administrator if problem persists"
        ]
        
        super().__init__(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_type=ErrorCategory.TASK_ERROR,
            message=f"Failed to submit {task_type} task",
            details=f"Could not submit task to processing queue: {error_details}",
            suggested_actions=suggested_actions,
            context={"task_type": task_type, "error_details": error_details}
        )


# ---------------------------------------------------------------------------
# Domain exception hierarchy (PR-8, Phase 3d)
#
# These are plain ``Exception`` subclasses — NOT ``HTTPException`` — so the
# service layer can raise domain-meaningful errors with no HTTP coupling. The
# global handler registered in ``app.main.create_app`` maps them back onto the
# exact same wire response FastAPI produces for ``HTTPException``:
#
#     JSONResponse(status_code=exc.status_code, content={"detail": exc.detail})
#
# so converting an ``HTTPException(status_code=S, detail=D)`` raise to the
# matching domain exception is response-preserving (byte-identical body + code).
# Each class carries a default ``status_code``; pass ``status_code=`` to cover
# the less common codes (413/502/503) while keeping the same semantic class.
# ---------------------------------------------------------------------------

class CellcraftError(Exception):
    """Base class for domain (non-HTTP) errors raised by the service layer."""

    status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR

    def __init__(self, detail: Any = None, *, status_code: Optional[int] = None):
        if status_code is not None:
            self.status_code = status_code
        self.detail = detail
        super().__init__(str(detail) if detail is not None else self.__class__.__name__)


class NotFoundError(CellcraftError):
    """Requested resource does not exist (maps to HTTP 404)."""

    status_code = status.HTTP_404_NOT_FOUND


class PermissionDeniedError(CellcraftError):
    """Caller is not allowed to perform the action (maps to HTTP 403)."""

    status_code = status.HTTP_403_FORBIDDEN


class ValidationFailedError(CellcraftError):
    """Input failed a business-rule validation (maps to HTTP 400)."""

    status_code = status.HTTP_400_BAD_REQUEST


class ExternalServiceError(CellcraftError):
    """A downstream service failed — Docker / GitHub / Redis (maps to HTTP 502)."""

    status_code = status.HTTP_502_BAD_GATEWAY


def log_error(error: Exception, context: Dict[str, Any] = None):
    """Log error with context information."""
    context_str = f" Context: {context}" if context else ""
    logger.error(f"Error occurred: {str(error)}{context_str}", exc_info=True)


def create_error_response(error: Exception, default_message: str = "An error occurred") -> JSONResponse:
    """Create a standardized error response from any exception."""
    if isinstance(error, CellCraftHTTPException):
        return JSONResponse(
            status_code=error.status_code,
            content=error.detail
        )
    elif isinstance(error, HTTPException):
        # Convert standard HTTPException to our format
        error_response = ErrorResponse(
            error_type=ErrorCategory.SERVER_ERROR,
            message=str(error.detail),
            details="An HTTP error occurred",
            suggested_actions=["Try again", "Contact support if problem persists"]
        )
        return JSONResponse(
            status_code=error.status_code,
            content=error_response.to_dict()
        )
    else:
        # Handle unexpected errors
        log_error(error)
        error_response = ErrorResponse(
            error_type=ErrorCategory.SERVER_ERROR,
            message=default_message,
            details="An unexpected error occurred on the server",
            suggested_actions=["Try again", "Contact support if problem persists"],
            context={"error_type": type(error).__name__}
        )
        return JSONResponse(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            content=error_response.to_dict()
        )