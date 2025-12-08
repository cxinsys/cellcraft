"""
Log archive utilities for task execution logs.

Provides functionality to archive execution logs per task_id,
enabling log preservation across multiple executions of the same algorithm node.
"""

import os
import shutil
from pathlib import Path
from typing import Dict, Optional, Any
from datetime import datetime


def archive_execution_logs(
    snakefile_dir: Path,
    task_id: str,
    include_meta: bool = True
) -> Dict[str, Any]:
    """
    Copy logs from current execution to executions/{task_id}/logs/.

    Archives the current execution logs to a task-specific directory,
    allowing previous execution logs to be preserved when the same
    algorithm node is re-executed with different settings.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        task_id: Unique task identifier for this execution
        include_meta: Whether to also copy meta.yml (default: True)

    Returns:
        Dict containing:
            - success: bool indicating if archival was successful
            - archived_path: Path to archived logs (if successful)
            - files_copied: List of files that were copied
            - error: Error message (if failed)

    Examples:
        >>> result = archive_execution_logs(Path("./user/john/workflow_1/algorithm_1"), "task-abc123")
        >>> result['success']
        True
        >>> result['archived_path']
        './user/john/workflow_1/algorithm_1/executions/task-abc123/logs'
    """
    result = {
        "success": False,
        "archived_path": None,
        "files_copied": [],
        "error": None,
        "timestamp": datetime.now().isoformat()
    }

    try:
        snakefile_dir = Path(snakefile_dir)

        # Source logs directory
        source_logs_dir = snakefile_dir / "logs"

        if not source_logs_dir.exists():
            result["error"] = f"Source logs directory does not exist: {source_logs_dir}"
            return result

        # Create executions/{task_id}/logs directory
        executions_dir = snakefile_dir / "executions" / task_id
        archived_logs_dir = executions_dir / "logs"
        archived_logs_dir.mkdir(parents=True, exist_ok=True, mode=0o777)

        # Copy all log files
        files_copied = []
        for log_file in source_logs_dir.iterdir():
            if log_file.is_file():
                dest_file = archived_logs_dir / log_file.name
                shutil.copy2(log_file, dest_file)
                files_copied.append(log_file.name)

        # Optionally copy meta.yml
        if include_meta:
            meta_file = snakefile_dir / "meta.yml"
            if meta_file.exists():
                dest_meta = executions_dir / "meta.yml"
                shutil.copy2(meta_file, dest_meta)
                files_copied.append("meta.yml (to execution root)")

        result["success"] = True
        result["archived_path"] = str(archived_logs_dir)
        result["files_copied"] = files_copied

        print(f"Archived {len(files_copied)} files to {archived_logs_dir}")

    except Exception as e:
        result["error"] = str(e)
        print(f"Failed to archive execution logs: {e}")

    return result


def get_execution_logs_path(
    base_path: str,
    username: str,
    workflow_id: int,
    algorithm_id: str,
    task_type: str,
    task_id: Optional[str] = None
) -> str:
    """
    Get logs path - archived if task_id provided and exists, else current.

    This function provides backward compatibility by checking for archived
    logs first (if task_id is provided) and falling back to the current
    logs directory if no archive exists.

    Args:
        base_path: Base path for user directories
        username: User's username
        workflow_id: Workflow database ID
        algorithm_id: Algorithm or visualization ID
        task_type: Task type ('visualization' or other)
        task_id: Optional task ID to look for archived logs

    Returns:
        str: Path to logs directory (archived if found, else current)

    Examples:
        >>> get_execution_logs_path("./user", "john", 1, "algo_1", "compile", "task-123")
        './user/john/workflow_1/algorithm_algo_1/executions/task-123/logs'
        >>> get_execution_logs_path("./user", "john", 1, "algo_1", "compile")
        './user/john/workflow_1/algorithm_algo_1/logs'
    """
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

    # Fallback to original path (backward compatibility)
    return f"{algo_base}/logs"


def find_archived_execution(
    snakefile_dir: Path,
    task_id: str
) -> Optional[Path]:
    """
    Find archived logs for a specific task_id.

    Checks if an archived execution exists for the given task_id
    within the snakefile directory.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        task_id: Task ID to look for

    Returns:
        Optional[Path]: Path to archived execution directory if found, else None

    Examples:
        >>> find_archived_execution(Path("./user/john/workflow_1/algorithm_1"), "task-123")
        PosixPath('./user/john/workflow_1/algorithm_1/executions/task-123')
    """
    snakefile_dir = Path(snakefile_dir)
    archived_dir = snakefile_dir / "executions" / task_id

    if archived_dir.exists() and archived_dir.is_dir():
        return archived_dir

    return None


def list_archived_executions(snakefile_dir: Path) -> list:
    """
    List all archived executions for an algorithm directory.

    Returns a list of all task_ids that have archived logs
    in the executions directory.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)

    Returns:
        list: List of task_ids with archived logs

    Examples:
        >>> list_archived_executions(Path("./user/john/workflow_1/algorithm_1"))
        ['task-123', 'task-456', 'task-789']
    """
    snakefile_dir = Path(snakefile_dir)
    executions_dir = snakefile_dir / "executions"

    if not executions_dir.exists():
        return []

    return [
        d.name for d in executions_dir.iterdir()
        if d.is_dir() and (d / "logs").exists()
    ]


def get_archived_logs_dir(
    snakefile_dir: Path,
    task_id: str
) -> Optional[Path]:
    """
    Get the archived logs directory path for a specific task.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        task_id: Task ID to get logs for

    Returns:
        Optional[Path]: Path to archived logs directory if exists, else None
    """
    archived_execution = find_archived_execution(snakefile_dir, task_id)
    if archived_execution:
        logs_dir = archived_execution / "logs"
        if logs_dir.exists():
            return logs_dir
    return None
