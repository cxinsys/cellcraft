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


def setup_execution_logs_symlink(
    snakefile_dir: Path,
    task_id: str
) -> Dict[str, Any]:
    """
    Set up logs/ as a symlink to executions/{task_id}/logs/ at task start.

    This function creates the execution-specific logs directory and sets up
    a symbolic link so that logs written to logs/ are actually stored in
    executions/{task_id}/logs/. This ensures logs are preserved even if
    the task is cancelled.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        task_id: Unique task identifier for this execution

    Returns:
        Dict containing:
            - success: bool indicating if setup was successful
            - logs_path: Path to actual logs directory (executions/{task_id}/logs/)
            - symlink_path: Path to symlink (logs/)
            - error: Error message (if failed)

    Examples:
        >>> result = setup_execution_logs_symlink(Path("./user/john/workflow_1/algorithm_1"), "task-abc123")
        >>> result['success']
        True
        >>> result['logs_path']
        './user/john/workflow_1/algorithm_1/executions/task-abc123/logs'
    """
    result = {
        "success": False,
        "logs_path": None,
        "symlink_path": None,
        "error": None,
        "timestamp": datetime.now().isoformat()
    }

    try:
        snakefile_dir = Path(snakefile_dir)

        # Create executions/{task_id}/logs directory (actual storage location)
        executions_dir = snakefile_dir / "executions" / task_id
        actual_logs_dir = executions_dir / "logs"
        actual_logs_dir.mkdir(parents=True, exist_ok=True, mode=0o777)

        # Path for the symlink
        symlink_path = snakefile_dir / "logs"

        # Remove existing logs/ if it exists (symlink or directory)
        if symlink_path.exists() or symlink_path.is_symlink():
            if symlink_path.is_symlink():
                # Remove existing symlink
                symlink_path.unlink()
                print(f"Removed existing symlink: {symlink_path}")
            elif symlink_path.is_dir():
                # Remove existing directory (previous execution logs already archived)
                shutil.rmtree(symlink_path)
                print(f"Removed existing logs directory: {symlink_path}")
            else:
                # Remove file if somehow it's a file
                symlink_path.unlink()
                print(f"Removed existing logs file: {symlink_path}")

        # Create symlink: logs/ -> executions/{task_id}/logs/
        # Use relative path for portability
        rel_target = os.path.relpath(actual_logs_dir, symlink_path.parent)
        symlink_path.symlink_to(rel_target)

        result["success"] = True
        result["logs_path"] = str(actual_logs_dir)
        result["symlink_path"] = str(symlink_path)

        print(f"Created logs symlink: {symlink_path} -> {rel_target}")

    except Exception as e:
        result["error"] = str(e)
        print(f"Failed to setup execution logs symlink: {e}")

    return result


def cleanup_task_results(
    snakefile_dir: Path,
    preserve_folder: bool = True
) -> Dict[str, Any]:
    """
    Clean up results folder contents when a task is cancelled.

    Removes all files inside the results/ folder while optionally
    preserving the folder structure. Handles symlinks carefully to
    preserve cache originals.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        preserve_folder: If True, keep the results/ folder but delete contents.
                        If False, delete the entire results/ folder.

    Returns:
        Dict containing:
            - success: bool indicating if cleanup was successful
            - files_removed: List of files that were removed
            - symlinks_removed: List of symlinks that were removed
            - error: Error message (if failed)

    Examples:
        >>> result = cleanup_task_results(Path("./user/john/workflow_1/algorithm_1"))
        >>> result['success']
        True
        >>> result['files_removed']
        ['output.h5ad', 'network.sif']
    """
    result = {
        "success": False,
        "files_removed": [],
        "symlinks_removed": [],
        "error": None,
        "timestamp": datetime.now().isoformat()
    }

    try:
        snakefile_dir = Path(snakefile_dir)
        results_dir = snakefile_dir / "results"

        if not results_dir.exists():
            result["success"] = True
            result["error"] = "Results directory does not exist (nothing to clean)"
            return result

        files_removed = []
        symlinks_removed = []

        # Iterate through all items in results directory
        for item in results_dir.iterdir():
            try:
                if item.is_symlink():
                    # Remove symlink only (preserve cache original)
                    item.unlink()
                    symlinks_removed.append(item.name)
                    print(f"Removed symlink: {item}")
                elif item.is_file():
                    # Remove regular file
                    item.unlink()
                    files_removed.append(item.name)
                    print(f"Removed file: {item}")
                elif item.is_dir():
                    # Remove subdirectory
                    shutil.rmtree(item)
                    files_removed.append(f"{item.name}/ (directory)")
                    print(f"Removed directory: {item}")
            except Exception as item_error:
                print(f"Warning: Failed to remove {item}: {item_error}")

        # Optionally remove the results folder itself
        if not preserve_folder and results_dir.exists():
            results_dir.rmdir()
            print(f"Removed results directory: {results_dir}")

        result["success"] = True
        result["files_removed"] = files_removed
        result["symlinks_removed"] = symlinks_removed

        total_removed = len(files_removed) + len(symlinks_removed)
        print(f"Cleaned up {total_removed} items from results directory")

    except Exception as e:
        result["error"] = str(e)
        print(f"Failed to cleanup task results: {e}")

    return result


def copy_meta_to_execution(
    snakefile_dir: Path,
    task_id: str
) -> Dict[str, Any]:
    """
    Copy meta.yml to executions/{task_id}/ directory.

    This is called after task completion to archive the meta.yml file
    alongside the logs in the execution-specific directory.

    Args:
        snakefile_dir: Directory containing the Snakefile (algorithm directory)
        task_id: Unique task identifier for this execution

    Returns:
        Dict containing:
            - success: bool indicating if copy was successful
            - dest_path: Path to copied meta.yml
            - error: Error message (if failed)
    """
    result = {
        "success": False,
        "dest_path": None,
        "error": None
    }

    try:
        snakefile_dir = Path(snakefile_dir)
        meta_file = snakefile_dir / "meta.yml"

        if not meta_file.exists():
            result["error"] = "meta.yml does not exist"
            return result

        executions_dir = snakefile_dir / "executions" / task_id
        if not executions_dir.exists():
            executions_dir.mkdir(parents=True, exist_ok=True, mode=0o777)

        dest_meta = executions_dir / "meta.yml"
        shutil.copy2(meta_file, dest_meta)

        result["success"] = True
        result["dest_path"] = str(dest_meta)
        print(f"Copied meta.yml to {dest_meta}")

    except Exception as e:
        result["error"] = str(e)
        print(f"Failed to copy meta.yml: {e}")

    return result
