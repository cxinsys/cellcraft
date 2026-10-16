"""
Path utility functions for secure file system operations.

Provides centralized path construction and validation to prevent
path traversal attacks and ensure consistent file access patterns.
"""

import os
from pathlib import Path
from typing import Optional


def get_logs_base_path() -> str:
    """
    Get base path for logs directory.

    Returns the base path for user logs. Can be overridden via
    CELLCRAFT_USER_BASE_PATH environment variable for testing.

    Returns:
        str: Base path for user directories (default: "./user")

    Examples:
        >>> get_logs_base_path()
        './user'
        >>> os.environ['CELLCRAFT_USER_BASE_PATH'] = '/tmp/test'
        >>> get_logs_base_path()
        '/tmp/test'
    """
    return os.getenv("CELLCRAFT_USER_BASE_PATH", "./user")


def construct_logs_path(
    username: str,
    workflow_id: int,
    algorithm_id: str,
    task_type: str,
    base_path: Optional[str] = None,
    task_id: Optional[str] = None
) -> str:
    """
    Construct logs directory path for a task.

    Creates a standardized path to task logs based on user, workflow,
    and task identifiers. Supports different task types (compile, visualization).

    If task_id is provided, checks for archived logs in executions/{task_id}/logs
    first, falling back to the current logs directory for backward compatibility.

    For visualization tasks, includes fallback logic to search other visualization_*
    folders when the primary path doesn't exist. This handles legacy tasks where
    algorithm_id was incorrectly stored (pointing to connected algorithm node
    instead of the visualization node itself).

    Args:
        username: User's username
        workflow_id: Workflow database ID
        algorithm_id: Algorithm or visualization ID
        task_type: Task type ('visualization' or other)
        base_path: Optional base path override (for testing)
        task_id: Optional task ID to look for archived logs

    Returns:
        str: Path to logs directory

    Examples:
        >>> construct_logs_path("user1", 123, "algo_1", "compile")
        './user/user1/workflow_123/algorithm_algo_1/logs'
        >>> construct_logs_path("user1", 123, "viz_1", "visualization")
        './user/user1/workflow_123/visualization_viz_1/logs'
        >>> construct_logs_path("user1", 123, "algo_1", "compile", task_id="task-abc")
        './user/user1/workflow_123/algorithm_algo_1/executions/task-abc/logs'

    Security:
        Does NOT validate path safety. Use is_safe_path() to verify.
    """
    if base_path is None:
        base_path = get_logs_base_path()

    # Determine base algorithm/visualization directory
    if task_type == 'visualization':
        algo_base = f"{base_path}/{username}/workflow_{workflow_id}/visualization_{algorithm_id}"
    else:
        algo_base = f"{base_path}/{username}/workflow_{workflow_id}/algorithm_{algorithm_id}"

    # Check archived path first if task_id provided
    if task_id:
        archived_path = f"{algo_base}/executions/{task_id}/logs"
        if os.path.exists(archived_path):
            return archived_path

    # Check default logs path
    default_path = f"{algo_base}/logs"
    if os.path.exists(default_path):
        return default_path

    # Fallback for legacy visualization tasks:
    # Search other visualization_* folders for archived logs with matching task_id
    # This handles cases where algorithm_id was incorrectly stored in legacy tasks
    if task_type == 'visualization' and task_id:
        found_path = _search_visualization_logs_fallback(
            base_path, username, workflow_id, task_id
        )
        if found_path:
            return found_path

    # Return default path even if it doesn't exist (caller handles 404)
    return default_path


def _search_visualization_logs_fallback(
    base_path: str,
    username: str,
    workflow_id: int,
    task_id: str
) -> Optional[str]:
    """
    Search for visualization logs in alternative folders.

    This fallback handles legacy visualization tasks where algorithm_id
    was incorrectly stored (pointing to connected algorithm node instead
    of the visualization node itself).

    Args:
        base_path: Base path for user directories
        username: User's username
        workflow_id: Workflow database ID
        task_id: Task ID to search for

    Returns:
        Optional[str]: Path to logs if found, None otherwise
    """
    workflow_base = f"{base_path}/{username}/workflow_{workflow_id}"

    if not os.path.exists(workflow_base):
        return None

    try:
        for folder in os.listdir(workflow_base):
            if folder.startswith('visualization_'):
                # Check archived logs path
                candidate_archived = f"{workflow_base}/{folder}/executions/{task_id}/logs"
                if os.path.exists(candidate_archived):
                    return candidate_archived

                # Check if current logs folder has run.log with matching task_id
                candidate_current = f"{workflow_base}/{folder}/logs"
                if os.path.exists(candidate_current):
                    run_log_path = f"{candidate_current}/run.log"
                    if os.path.exists(run_log_path):
                        try:
                            with open(run_log_path, 'r', encoding='utf-8') as f:
                                # Read first few KB to check for task_id
                                content = f.read(4096)
                                if task_id in content:
                                    return candidate_current
                        except (IOError, OSError):
                            # Skip files that can't be read
                            pass
    except (IOError, OSError):
        # Directory listing failed
        pass

    return None


def is_safe_path(path: str, base_dir: str) -> bool:
    """
    Validate that a path is safe (no path traversal).

    Ensures that the resolved absolute path is within the base directory,
    preventing path traversal attacks using ../ or symlinks.

    Args:
        path: Path to validate
        base_dir: Base directory that path must be within

    Returns:
        bool: True if path is safe, False if path traversal detected

    Examples:
        >>> is_safe_path("./user/john/logs", "./user")
        True
        >>> is_safe_path("./user/../../../etc/passwd", "./user")
        False
        >>> is_safe_path("/tmp/symlink_to_etc", "./user")
        False

    Security:
        - Resolves symlinks to detect symlink attacks
        - Converts to absolute paths for comparison
        - Checks if path is within base_dir tree
    """
    try:
        # Resolve to absolute paths, following symlinks
        abs_path = Path(path).resolve()
        abs_base = Path(base_dir).resolve()

        # Check if path is within base directory
        return abs_base in abs_path.parents or abs_path == abs_base
    except (ValueError, RuntimeError):
        # Path resolution failed (e.g., broken symlink, invalid path)
        return False


def sanitize_filename(filename: str) -> Optional[str]:
    """
    Sanitize a filename to prevent path traversal.

    Removes directory separators and path traversal sequences
    from filename strings.

    Args:
        filename: Filename to sanitize

    Returns:
        Optional[str]: Sanitized filename, or None if invalid

    Examples:
        >>> sanitize_filename("safe.log")
        'safe.log'
        >>> sanitize_filename("../../../etc/passwd")
        None
        >>> sanitize_filename("..\\..\\windows\\system.ini")
        None

    Security:
        Rejects filenames containing:
        - Directory separators (/ or \\)
        - Parent directory references (..)
        - Null bytes
    """
    if not filename:
        return None

    # Reject path traversal attempts
    dangerous_patterns = ['..', '/', '\\', '\0']
    for pattern in dangerous_patterns:
        if pattern in filename:
            return None

    # Additional validation: must be valid filename
    if not filename.strip() or filename.startswith('.'):
        return None

    return filename


def validate_export_filename(filename: str, task_id: str) -> bool:
    """
    Validate filename for log export endpoints.

    Ensures filename is safe for use in export operations.
    Used by /export/txt/{filename} and /export/json/{filename} endpoints.

    Args:
        filename: Filename to validate
        task_id: Associated task ID (for logging)

    Returns:
        bool: True if filename is safe for export

    Examples:
        >>> validate_export_filename("run.log", "task-123")
        True
        >>> validate_export_filename("../../../etc/passwd", "task-123")
        False

    Security:
        Uses sanitize_filename() to prevent path traversal.
    """
    sanitized = sanitize_filename(filename)
    return sanitized is not None and sanitized == filename
