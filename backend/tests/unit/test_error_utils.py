"""
Unit tests for error handling utilities.

Test coverage:
- ErrorCategory enum
- ErrorResponse class
- CellCraftHTTPException base class
- Specialized exception classes
- Utility functions (log_error, create_error_response)
"""
import pytest
import json
from unittest.mock import patch
from fastapi import HTTPException, status
from fastapi.responses import JSONResponse

from app.common.utils.error_utils import (
    ErrorCategory,
    ErrorResponse,
    CellCraftHTTPException,
    PluginNotFoundError,
    ScriptNotFoundError,
    FileNotFoundError,
    ValidationError,
    WorkflowError,
    SnakefileGenerationError,
    TaskSubmissionError,
    log_error,
    create_error_response
)


@pytest.mark.unit
@pytest.mark.error
class TestErrorCategory:
    """Test ErrorCategory enum."""

    def test_error_category_values(self):
        """Test all ErrorCategory enum values."""
        # Assert all 8 error categories exist with correct values
        assert ErrorCategory.VALIDATION == "validation"
        assert ErrorCategory.NOT_FOUND == "not_found"
        assert ErrorCategory.PERMISSION == "permission"
        assert ErrorCategory.SERVER_ERROR == "server_error"
        assert ErrorCategory.PLUGIN_ERROR == "plugin_error"
        assert ErrorCategory.FILE_ERROR == "file_error"
        assert ErrorCategory.WORKFLOW_ERROR == "workflow_error"
        assert ErrorCategory.TASK_ERROR == "task_error"

        # Verify we have exactly 8 categories
        assert len(ErrorCategory) == 8


@pytest.mark.unit
@pytest.mark.error
class TestErrorResponse:
    """Test ErrorResponse class."""

    def test_error_response_initialization(self):
        """Test ErrorResponse initialization with all fields."""
        # Arrange
        error_type = ErrorCategory.VALIDATION
        message = "Test error message"
        details = "Detailed error information"
        suggested_actions = ["Action 1", "Action 2"]
        context = {"field": "test_field", "value": "test_value"}

        # Act
        error_response = ErrorResponse(
            error_type=error_type,
            message=message,
            details=details,
            suggested_actions=suggested_actions,
            context=context
        )

        # Assert
        assert error_response.error_type == error_type
        assert error_response.message == message
        assert error_response.details == details
        assert error_response.suggested_actions == suggested_actions
        assert error_response.context == context

    def test_error_response_to_dict(self):
        """Test ErrorResponse to_dict method output format."""
        # Arrange
        error_response = ErrorResponse(
            error_type=ErrorCategory.NOT_FOUND,
            message="Resource not found",
            details="The requested resource does not exist",
            suggested_actions=["Check the resource ID", "Try again"],
            context={"resource_id": "123"}
        )

        # Act
        result = error_response.to_dict()

        # Assert
        assert "error" in result
        assert result["error"]["type"] == "not_found"
        assert result["error"]["message"] == "Resource not found"
        assert result["error"]["details"] == "The requested resource does not exist"
        assert len(result["error"]["suggested_actions"]) == 2
        assert result["error"]["context"]["resource_id"] == "123"

    def test_error_response_with_minimal_parameters(self):
        """Test ErrorResponse with only required parameters."""
        # Arrange & Act
        error_response = ErrorResponse(
            error_type=ErrorCategory.VALIDATION,
            message="Validation failed"
        )

        # Assert
        assert error_response.error_type == ErrorCategory.VALIDATION
        assert error_response.message == "Validation failed"
        assert error_response.details is None
        assert error_response.suggested_actions == []
        assert error_response.context == {}

    def test_error_response_to_dict_with_empty_optionals(self):
        """Test to_dict handles empty optional fields correctly."""
        # Arrange
        error_response = ErrorResponse(
            error_type=ErrorCategory.SERVER_ERROR,
            message="Internal error"
        )

        # Act
        result = error_response.to_dict()

        # Assert
        assert "error" in result
        assert result["error"]["type"] == "server_error"
        assert result["error"]["message"] == "Internal error"
        assert result["error"]["details"] is None
        assert result["error"]["suggested_actions"] == []
        assert result["error"]["context"] == {}

    def test_error_response_default_initialization_with_none_values(self):
        """Test ErrorResponse converts None optional parameters to proper defaults.

        When optional parameters are explicitly set to None:
        - suggested_actions should become [] (empty list, not None)
        - context should become {} (empty dict, not None)

        This prevents NoneType errors in consumers expecting list/dict.
        """
        # Arrange & Act - Explicitly pass None for optional parameters
        error_response_1 = ErrorResponse(
            error_type=ErrorCategory.VALIDATION,
            message="Test error 1",
            details="Test details",
            suggested_actions=None,
            context=None
        )

        # Assert - Verify None converts to proper defaults
        assert error_response_1.suggested_actions == []
        assert error_response_1.suggested_actions is not None
        assert isinstance(error_response_1.suggested_actions, list)

        assert error_response_1.context == {}
        assert error_response_1.context is not None
        assert isinstance(error_response_1.context, dict)

        # Arrange & Act - Create second instance to ensure no shared mutable defaults
        error_response_2 = ErrorResponse(
            error_type=ErrorCategory.NOT_FOUND,
            message="Test error 2",
            suggested_actions=None,
            context=None
        )

        # Act - Modify second instance
        error_response_2.suggested_actions.append("Action for instance 2")
        error_response_2.context["key"] = "value for instance 2"

        # Assert - Verify first instance is not affected (no shared mutable defaults)
        assert len(error_response_1.suggested_actions) == 0
        assert len(error_response_1.context) == 0
        assert len(error_response_2.suggested_actions) == 1
        assert len(error_response_2.context) == 1


@pytest.mark.unit
@pytest.mark.error
class TestCellCraftHTTPException:
    """Test CellCraftHTTPException base class."""

    def test_cellcraft_http_exception_creates_error_response(self):
        """Test CellCraftHTTPException creates ErrorResponse from parameters."""
        # Arrange & Act
        exception = CellCraftHTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            error_type=ErrorCategory.VALIDATION,
            message="Custom error message",
            details="Custom details",
            suggested_actions=["Action 1"],
            context={"key": "value"}
        )

        # Assert - Verify error_response is created correctly
        assert exception.status_code == status.HTTP_400_BAD_REQUEST
        assert exception.error_response.error_type == ErrorCategory.VALIDATION
        assert exception.error_response.message == "Custom error message"
        assert exception.error_response.details == "Custom details"
        assert exception.error_response.suggested_actions == ["Action 1"]
        assert exception.error_response.context == {"key": "value"}

    def test_cellcraft_http_exception_inherits_from_http_exception(self):
        """Test CellCraftHTTPException inherits from HTTPException with correct status code."""
        # Arrange & Act
        exception = CellCraftHTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_type=ErrorCategory.SERVER_ERROR,
            message="Server error"
        )

        # Assert - Verify it's an HTTPException
        assert isinstance(exception, HTTPException)
        assert exception.status_code == status.HTTP_500_INTERNAL_SERVER_ERROR

    def test_cellcraft_http_exception_detail_contains_error_dict(self):
        """Test CellCraftHTTPException detail contains structured error dictionary."""
        # Arrange & Act
        exception = CellCraftHTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            error_type=ErrorCategory.PLUGIN_ERROR,
            message="Plugin initialization failed",
            details="Plugin 'test' could not be loaded",
            suggested_actions=["Check plugin configuration", "Verify dependencies"],
            context={"plugin_name": "test", "error_code": "INIT_FAILED"}
        )

        # Assert - Verify detail is the error response dictionary
        assert isinstance(exception.detail, dict)
        assert "error" in exception.detail
        assert exception.detail["error"]["type"] == "plugin_error"
        assert exception.detail["error"]["message"] == "Plugin initialization failed"
        assert exception.detail["error"]["context"]["plugin_name"] == "test"
        assert exception.detail["error"]["context"]["error_code"] == "INIT_FAILED"

    def test_cellcraft_http_exception_with_all_error_categories(self):
        """Test CellCraftHTTPException works correctly with all 8 ErrorCategory types.

        Verifies that:
        - Each ErrorCategory can be used to create an exception
        - error_type propagates correctly to error_response
        - detail dict structure remains consistent across all categories
        """
        # Arrange - Get all 8 ErrorCategory values
        all_categories = [
            ErrorCategory.VALIDATION,
            ErrorCategory.NOT_FOUND,
            ErrorCategory.PERMISSION,
            ErrorCategory.SERVER_ERROR,
            ErrorCategory.PLUGIN_ERROR,
            ErrorCategory.FILE_ERROR,
            ErrorCategory.WORKFLOW_ERROR,
            ErrorCategory.TASK_ERROR
        ]

        # Act & Assert - Loop through all categories
        for category in all_categories:
            # Act - Create exception with current category
            exception = CellCraftHTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                error_type=category,
                message=f"Test message for {category.value}",
                details="Test details",
                suggested_actions=["Action 1", "Action 2"],
                context={"test_key": "test_value"}
            )

            # Assert - error_type propagates correctly
            assert exception.error_response.error_type == category
            assert exception.error_response.error_type.value == category.value

            # Assert - detail dict has consistent 'error' key structure
            assert isinstance(exception.detail, dict)
            assert "error" in exception.detail
            assert exception.detail["error"]["type"] == category.value
            assert exception.detail["error"]["message"] == f"Test message for {category.value}"
            assert exception.detail["error"]["details"] == "Test details"
            assert exception.detail["error"]["suggested_actions"] == ["Action 1", "Action 2"]
            assert exception.detail["error"]["context"] == {"test_key": "test_value"}


@pytest.mark.unit
@pytest.mark.error
class TestPluginNotFoundError:
    """Test PluginNotFoundError exception."""

    def test_plugin_not_found_error_basic(self):
        """Test PluginNotFoundError with plugin name only."""
        # Arrange
        plugin_name = "test_plugin"

        # Act
        error = PluginNotFoundError(plugin_name)

        # Assert
        assert error.status_code == status.HTTP_404_NOT_FOUND
        assert error.error_response.error_type == ErrorCategory.NOT_FOUND
        assert plugin_name in error.error_response.message
        assert len(error.error_response.suggested_actions) == 3
        assert error.error_response.context["plugin_name"] == plugin_name

    def test_plugin_not_found_error_with_available_plugins(self):
        """Test PluginNotFoundError with available plugins list."""
        # Arrange
        plugin_name = "missing_plugin"
        available_plugins = ["plugin1", "plugin2", "plugin3"]

        # Act
        error = PluginNotFoundError(plugin_name, available_plugins)

        # Assert
        assert error.error_response.context["available_plugins"] == available_plugins
        # Should have 4 suggested actions (3 default + 1 with available plugins)
        assert len(error.error_response.suggested_actions) == 4


@pytest.mark.unit
@pytest.mark.error
class TestScriptNotFoundError:
    """Test ScriptNotFoundError exception."""

    def test_script_not_found_error(self):
        """Test ScriptNotFoundError with script and plugin names."""
        # Arrange
        script_name = "test_script.py"
        plugin_name = "test_plugin"
        available_scripts = ["script1.py", "script2.py"]

        # Act
        error = ScriptNotFoundError(script_name, plugin_name, available_scripts)

        # Assert
        assert error.status_code == status.HTTP_404_NOT_FOUND
        assert error.error_response.error_type == ErrorCategory.NOT_FOUND
        assert script_name in error.error_response.message
        assert plugin_name in error.error_response.message
        assert error.error_response.context["script_name"] == script_name
        assert error.error_response.context["plugin_name"] == plugin_name
        assert error.error_response.context["available_scripts"] == available_scripts


@pytest.mark.unit
@pytest.mark.error
class TestFileNotFoundError:
    """Test FileNotFoundError exception."""

    def test_file_not_found_error_default_type(self):
        """Test FileNotFoundError with default file type."""
        # Arrange
        file_path = "/path/to/missing/file.txt"

        # Act
        error = FileNotFoundError(file_path)

        # Assert
        assert error.status_code == status.HTTP_404_NOT_FOUND
        assert error.error_response.error_type == ErrorCategory.FILE_ERROR
        assert file_path in error.error_response.details
        assert error.error_response.context["file_path"] == file_path
        assert error.error_response.context["file_type"] == "file"

    def test_file_not_found_error_custom_type(self):
        """Test FileNotFoundError with custom file type."""
        # Arrange
        file_path = "/data/dataset.h5ad"
        file_type = "dataset"

        # Act
        error = FileNotFoundError(file_path, file_type)

        # Assert
        assert error.error_response.context["file_type"] == file_type
        assert "Dataset" in error.error_response.message  # Title case


@pytest.mark.unit
@pytest.mark.error
class TestValidationError:
    """Test ValidationError exception."""

    def test_validation_error_basic(self):
        """Test ValidationError with basic parameters."""
        # Arrange
        field_name = "email"
        message = "Invalid email format"

        # Act
        error = ValidationError(field_name, message)

        # Assert
        assert error.status_code == status.HTTP_400_BAD_REQUEST
        assert error.error_response.error_type == ErrorCategory.VALIDATION
        assert field_name in error.error_response.message
        assert error.error_response.details == message
        assert error.error_response.context["field_name"] == field_name

    def test_validation_error_with_required_format(self):
        """Test ValidationError with required format specification."""
        # Arrange
        field_name = "date"
        message = "Invalid date format"
        required_format = "YYYY-MM-DD"

        # Act
        error = ValidationError(field_name, message, required_format)

        # Assert
        assert error.error_response.context["required_format"] == required_format
        # Should have extra suggested action with required format
        assert any(required_format in action for action in error.error_response.suggested_actions)


@pytest.mark.unit
@pytest.mark.error
class TestWorkflowError:
    """Test WorkflowError exception."""

    def test_workflow_error(self):
        """Test WorkflowError with workflow ID and message."""
        # Arrange
        workflow_id = "workflow_123"
        message = "Workflow compilation failed"
        error_type = "compilation"

        # Act
        error = WorkflowError(workflow_id, message, error_type)

        # Assert
        assert error.status_code == status.HTTP_400_BAD_REQUEST
        assert error.error_response.error_type == ErrorCategory.WORKFLOW_ERROR
        assert error_type in error.error_response.message
        assert error.error_response.details == message
        assert error.error_response.context["workflow_id"] == workflow_id
        assert error.error_response.context["error_type"] == error_type


@pytest.mark.unit
@pytest.mark.error
class TestSnakefileGenerationError:
    """Test SnakefileGenerationError exception."""

    def test_snakefile_generation_error(self):
        """Test SnakefileGenerationError with plugin and script details."""
        # Arrange
        plugin_name = "test_plugin"
        script_name = "visualization.py"
        error_details = "Missing required parameter 'input_file'"

        # Act
        error = SnakefileGenerationError(plugin_name, script_name, error_details)

        # Assert
        assert error.status_code == status.HTTP_500_INTERNAL_SERVER_ERROR
        assert error.error_response.error_type == ErrorCategory.SERVER_ERROR
        assert "Snakefile" in error.error_response.message
        assert plugin_name in error.error_response.details
        assert script_name in error.error_response.details
        assert error_details in error.error_response.details
        assert error.error_response.context["plugin_name"] == plugin_name
        assert error.error_response.context["script_name"] == script_name


@pytest.mark.unit
@pytest.mark.error
class TestTaskSubmissionError:
    """Test TaskSubmissionError exception."""

    def test_task_submission_error(self):
        """Test TaskSubmissionError with task type and details."""
        # Arrange
        task_type = "workflow_execution"
        error_details = "Queue is full"

        # Act
        error = TaskSubmissionError(task_type, error_details)

        # Assert
        assert error.status_code == status.HTTP_500_INTERNAL_SERVER_ERROR
        assert error.error_response.error_type == ErrorCategory.TASK_ERROR
        assert task_type in error.error_response.message
        assert error_details in error.error_response.details
        assert error.error_response.context["task_type"] == task_type
        assert error.error_response.context["error_details"] == error_details


@pytest.mark.unit
@pytest.mark.error
class TestLogError:
    """Test log_error utility function."""

    @patch('app.common.utils.error_utils.logger')
    def test_log_error_with_context(self, mock_logger):
        """Test log_error function with context and exact kwargs validation."""
        # Arrange
        error = Exception("Test error")
        context = {"user_id": "123", "action": "upload"}

        # Act
        log_error(error, context)

        # Assert - Validate exact call arguments
        mock_logger.error.assert_called_once()
        call_args, call_kwargs = mock_logger.error.call_args

        # Verify message contains error and context
        log_message = call_args[0]
        assert "Test error" in log_message
        assert "Context:" in log_message
        assert "user_id" in log_message
        assert "123" in log_message
        assert "action" in log_message
        assert "upload" in log_message

        # Verify exc_info is True for stack trace
        assert call_kwargs["exc_info"] is True

    @patch('app.common.utils.error_utils.logger')
    def test_log_error_without_context(self, mock_logger):
        """Test log_error function without context."""
        # Arrange
        error = ValueError("Invalid value provided")

        # Act
        log_error(error)

        # Assert
        mock_logger.error.assert_called_once()
        call_args, call_kwargs = mock_logger.error.call_args

        # Verify message contains error
        log_message = call_args[0]
        assert "Invalid value provided" in log_message

        # Verify exc_info is True even without context
        assert call_kwargs["exc_info"] is True

    @patch('app.common.utils.error_utils.logger')
    def test_log_error_exc_info_always_true(self, mock_logger):
        """Test log_error always passes exc_info=True for stack trace capture.

        Stack traces are critical for debugging, so exc_info must ALWAYS be True
        regardless of whether context is provided or error type.

        Tests both scenarios:
        - With context dict
        - Without context (None/omitted)
        """
        # Arrange - Test with different exception types
        error_with_context = ValueError("Value error with context")
        error_without_context = RuntimeError("Runtime error without context")
        context = {"operation": "data_processing", "user_id": "12345"}

        # Act - Test with context
        log_error(error_with_context, context)

        # Assert - Verify exc_info=True with context
        assert mock_logger.error.call_count == 1
        call_args_with, call_kwargs_with = mock_logger.error.call_args
        assert call_kwargs_with["exc_info"] is True
        assert "Value error with context" in call_args_with[0]
        assert "Context:" in call_args_with[0]

        # Reset mock for next test
        mock_logger.reset_mock()

        # Act - Test without context
        log_error(error_without_context)

        # Assert - Verify exc_info=True without context
        assert mock_logger.error.call_count == 1
        call_args_without, call_kwargs_without = mock_logger.error.call_args
        assert call_kwargs_without["exc_info"] is True
        assert "Runtime error without context" in call_args_without[0]

        # Reset mock for verification with explicit None
        mock_logger.reset_mock()

        # Act - Test with explicit None context
        log_error(error_with_context, None)

        # Assert - Verify exc_info=True with explicit None
        assert mock_logger.error.call_count == 1
        call_args_none, call_kwargs_none = mock_logger.error.call_args
        assert call_kwargs_none["exc_info"] is True


@pytest.mark.unit
@pytest.mark.error
class TestCreateErrorResponse:
    """Test create_error_response utility function."""

    def test_create_error_response_with_cellcraft_exception(self):
        """Test create_error_response with CellCraftHTTPException validates JSON structure."""
        # Arrange
        error = PluginNotFoundError("test_plugin")

        # Act
        response = create_error_response(error)

        # Assert - Validate response type and status
        assert isinstance(response, JSONResponse)
        assert response.status_code == status.HTTP_404_NOT_FOUND

        # Deserialize and validate JSON structure
        response_data = json.loads(response.body.decode())
        assert "error" in response_data
        assert "type" in response_data["error"]
        assert "message" in response_data["error"]
        assert "details" in response_data["error"]
        assert "suggested_actions" in response_data["error"]
        assert "context" in response_data["error"]

        # Validate specific content
        assert response_data["error"]["type"] == ErrorCategory.NOT_FOUND
        assert "test_plugin" in response_data["error"]["message"]
        assert isinstance(response_data["error"]["suggested_actions"], list)
        assert "plugin_name" in response_data["error"]["context"]

    def test_create_error_response_json_schema_completeness(self):
        """Test create_error_response returns complete JSON schema with all required fields.

        Validates that the JSON response structure contains:
        - Top-level 'error' key
        - All 5 required nested fields: type, message, details, suggested_actions, context
        - Correct data types for each field (str, list, dict)
        """
        # Arrange
        error = PluginNotFoundError("test_plugin")

        # Act
        response = create_error_response(error)

        # Assert - Deserialize JSON body
        response_data = json.loads(response.body.decode())

        # Assert - Top-level 'error' key exists
        assert "error" in response_data
        error_obj = response_data["error"]

        # Assert - All 5 required fields exist
        assert "type" in error_obj
        assert "message" in error_obj
        assert "details" in error_obj
        assert "suggested_actions" in error_obj
        assert "context" in error_obj

        # Assert - Correct data types for each field
        assert isinstance(error_obj["type"], str)
        assert isinstance(error_obj["message"], str)
        # details can be str or None
        assert error_obj["details"] is None or isinstance(error_obj["details"], str)
        assert isinstance(error_obj["suggested_actions"], list)
        assert isinstance(error_obj["context"], dict)

    def test_create_error_response_generic_exception_includes_error_type(self):
        """Test create_error_response with generic Exception includes error type in context.

        Generic exceptions should:
        - Log the error via log_error()
        - Return 500 Internal Server Error status
        - Include exception class name in context['error_type']
        - Use default_message parameter if provided
        """
        # Arrange
        error = ValueError("Invalid input value")
        default_message = "Validation processing failed"

        # Act - Test with default_message
        with patch('app.common.utils.error_utils.log_error') as mock_log_error:
            response = create_error_response(error, default_message)

            # Assert - Verify log_error was called
            mock_log_error.assert_called_once_with(error)

        # Assert - Verify response status and type
        assert isinstance(response, JSONResponse)
        assert response.status_code == status.HTTP_500_INTERNAL_SERVER_ERROR

        # Assert - Verify context contains error_type
        response_data = json.loads(response.body.decode())
        assert "error" in response_data
        assert "context" in response_data["error"]
        assert "error_type" in response_data["error"]["context"]
        assert response_data["error"]["context"]["error_type"] == "ValueError"

        # Assert - Verify default_message is used
        assert response_data["error"]["message"] == default_message

        # Act - Test without default_message
        with patch('app.common.utils.error_utils.log_error') as mock_log_error:
            response_no_default = create_error_response(error)

            # Assert - Verify log_error was called
            mock_log_error.assert_called_once_with(error)

        # Assert - Verify default message is used when not provided
        response_no_default_data = json.loads(response_no_default.body.decode())
        assert response_no_default_data["error"]["message"] == "An error occurred"

    def test_create_error_response_with_http_exception(self):
        """Test create_error_response converts HTTPException to structured format."""
        # Arrange
        error = HTTPException(status_code=400, detail="Bad request")

        # Act
        response = create_error_response(error)

        # Assert - Validate response type and status
        assert isinstance(response, JSONResponse)
        assert response.status_code == 400

        # Deserialize and validate JSON structure - HTTPException is converted to structured format
        response_data = json.loads(response.body.decode())
        assert "error" in response_data
        assert "type" in response_data["error"]
        assert "message" in response_data["error"]
        assert "details" in response_data["error"]
        assert response_data["error"]["type"] == "server_error"
        assert response_data["error"]["message"] == "Bad request"
        assert response_data["error"]["details"] == "An HTTP error occurred"

    @patch('app.common.utils.error_utils.log_error')
    def test_create_error_response_with_generic_exception(self, mock_log_error):
        """Test create_error_response with generic Exception."""
        # Arrange
        error = ValueError("Invalid value")
        default_message = "An error occurred"

        # Act
        response = create_error_response(error, default_message)

        # Assert
        assert isinstance(response, JSONResponse)
        assert response.status_code == status.HTTP_500_INTERNAL_SERVER_ERROR
        mock_log_error.assert_called_once_with(error)
