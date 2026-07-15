"""
Admin domain service layer.

Extracted from ``admin/router.py`` in PR-8 (Phase 3d). This module holds the
business logic that previously lived inline in the endpoints: superuser-gated
listing/counting of users/files/workflows/tasks/plugins, user/file/workflow/task
mutations, the Docker-based system-stats aggregation, plugin synchronization
management, and container-manager administration. The router now only parses
requests, resolves dependencies, calls a service function, and returns its
result.

Exception policy (PR-8): the repeated ``is_superuser`` gate now raises
``PermissionDeniedError`` (-> 403) via ``_require_admin``, and the "not found"
guards raise ``NotFoundError`` (-> 404). The global handler in ``app.main`` maps
these onto the exact ``{"detail": ...}`` wire format, so responses are unchanged.

NOTE (domain boundaries): ``admin/crud.py`` reads other domains' tables directly
(users/files/workflows/tasks/plugins) for the admin dashboards. That cross-domain
data access is pre-existing and left as-is per the PR-8 scope note (boundary
hardening is out of scope here). The plugin-sync / container-manager 500 handlers
that swallow-and-stringify their inner errors are kept as ``HTTPException``
verbatim to guarantee byte-identical error bodies.
"""
import logging

from fastapi import HTTPException
from sqlalchemy.orm import Session

from app import models
from app.admin.schemas import Conditions
from app.admin import crud as crud_admin
from app.core.exceptions import (
    NotFoundError, PermissionDeniedError, ExternalServiceError
)
from app.plugin.sync import PluginSyncManager
from app.plugin.versioning import PluginVersionValidator
from app.shared.docker import container_manager

logger = logging.getLogger(__name__)


def _require_admin(current_user: models.User) -> None:
    """Raise 403 unless the caller is a superuser (admin-only gate)."""
    if not current_user.is_superuser:
        raise PermissionDeniedError("Access denied: Admins only")


def get_filtered_users(*, db: Session, current_user: models.User, amount: int,
                       page_num: int, sort: str, order: str, searchTerm: str):
    """Return a filtered/paginated user list (admin only; 404 if none match)."""
    _require_admin(current_user)

    conditions = Conditions(
        amount=amount, page_num=page_num, sort=sort, order=order, searchTerm=searchTerm
    )
    users, total_count = crud_admin.get_filtered_users(db, conditions)

    if not users:
        raise NotFoundError("Users not found")
    return users


def get_users_count(*, db: Session, current_user: models.User):
    """Return the total user count (admin only)."""
    _require_admin(current_user)
    return crud_admin.get_users_count(db)


def get_filtered_files(*, db: Session, current_user: models.User, amount: int,
                       page_num: int, sort: str, order: str, searchTerm: str):
    """Return a filtered/paginated file list with owner usernames (admin only)."""
    _require_admin(current_user)

    conditions = Conditions(
        amount=amount, page_num=page_num, sort=sort, order=order, searchTerm=searchTerm
    )
    files, total_count = crud_admin.get_filtered_files(db, conditions)

    if not files:
        raise NotFoundError("Files not found")

    formatted_files = []
    for file, username in files:
        formatted_file = {
            'id': file.id,
            'file_name': file.file_name,
            'file_path': file.file_path,
            'file_size': file.file_size,
            'folder': file.folder,
            'username': username,
            'user_id': file.user_id
        }
        formatted_files.append(formatted_file)

    return {
        'data': formatted_files,
        'total_count': total_count
    }


def get_files_count(*, db: Session, current_user: models.User):
    """Return the total file count (admin only)."""
    _require_admin(current_user)
    return crud_admin.get_files_count(db)


def get_filtered_workflows(*, db: Session, current_user: models.User, amount: int,
                           page_num: int, sort: str, order: str, searchTerm: str):
    """Return a filtered/paginated workflow list with owner usernames (admin only)."""
    _require_admin(current_user)

    conditions = Conditions(
        amount=amount, page_num=page_num, sort=sort, order=order, searchTerm=searchTerm
    )
    workflows, total_count = crud_admin.get_filtered_workflows(db, conditions)

    if not workflows:
        raise NotFoundError("Workflows not found")

    formatted_workflows = []
    for workflow, username in workflows:
        formatted_workflow = {
            'id': workflow.id,
            'title': workflow.title,
            'username': username,
            'updated_at': workflow.updated_at,
            'user_id': workflow.user_id
        }
        formatted_workflows.append(formatted_workflow)

    return {
        'data': formatted_workflows,
        'total_count': total_count
    }


def get_workflows_count(*, db: Session, current_user: models.User):
    """Return the total workflow count (admin only)."""
    _require_admin(current_user)
    return crud_admin.get_workflows_count(db)


def get_filtered_tasks(*, db: Session, current_user: models.User, amount: int,
                       page_num: int, sort: str, order: str, searchTerm: str):
    """Return a filtered/paginated task list with owner + workflow titles (admin only)."""
    _require_admin(current_user)

    conditions = Conditions(
        amount=amount, page_num=page_num, sort=sort, order=order, searchTerm=searchTerm
    )
    tasks, total_count = crud_admin.get_filtered_tasks(db, conditions)

    if not tasks:
        raise NotFoundError("Tasks not found")

    formatted_tasks = []
    for task, username, workflow_title in tasks:
        formatted_task = {
            'id': task.id,
            'user_id': task.user_id,
            'workflow_id': task.workflow_id,
            'username': username,
            'workflow_title': workflow_title,
            'status': task.status,
            'start_time': task.start_time
        }
        formatted_tasks.append(formatted_task)

    return {
        'data': formatted_tasks,
        'total_count': total_count
    }


def get_tasks_count(*, db: Session, current_user: models.User):
    """Return the total task count (admin only)."""
    _require_admin(current_user)
    return crud_admin.get_tasks_count(db)


def get_filtered_plugins(*, db: Session, current_user: models.User, amount: int,
                         page_num: int, sort: str, order: str, searchTerm: str):
    """Return a filtered/paginated plugin list (admin only; 404 if none match)."""
    _require_admin(current_user)

    conditions = Conditions(
        amount=amount, page_num=page_num, sort=sort, order=order, searchTerm=searchTerm
    )
    plugins, total_count = crud_admin.get_filtered_plugins(db, conditions)

    if not plugins:
        raise NotFoundError("Plugins not found")
    return plugins


def get_plugins_count(*, db: Session, current_user: models.User):
    """Return the total plugin count (admin only)."""
    _require_admin(current_user)
    return crud_admin.get_plugins_count(db)


def update_user(*, db: Session, current_user: models.User, user_id: int, user_data: dict):
    """Update a user record (admin only; 404 if the user is missing)."""
    _require_admin(current_user)

    user = crud_admin.update_user(db, user_id, user_data)
    if not user:
        raise NotFoundError("User not found")
    return user


def delete_user(*, db: Session, current_user: models.User, user_id: int) -> dict:
    """Delete a user (admin only; 404 if the user is missing)."""
    _require_admin(current_user)

    success = crud_admin.delete_user(db, user_id)
    if not success:
        raise NotFoundError("User not found")
    return {"message": "User deleted successfully"}


def delete_file(*, db: Session, current_user: models.User, file_id: int) -> dict:
    """Delete a file record (admin only; 404 if the file is missing)."""
    _require_admin(current_user)

    success = crud_admin.delete_file(db, file_id)
    if not success:
        raise NotFoundError("File not found")
    return {"message": "File deleted successfully"}


def delete_workflow(*, db: Session, current_user: models.User, workflow_id: int) -> dict:
    """Delete a workflow record (admin only; 404 if the workflow is missing)."""
    _require_admin(current_user)

    success = crud_admin.delete_workflow(db, workflow_id)
    if not success:
        raise NotFoundError("Workflow not found")
    return {"message": "Workflow deleted successfully"}


def cancel_task(*, db: Session, current_user: models.User, task_id: int) -> dict:
    """Cancel a task (admin only; 404 if the task is missing)."""
    _require_admin(current_user)

    success = crud_admin.cancel_task(db, task_id)
    if not success:
        raise NotFoundError("Task not found")
    return {"message": "Task cancelled successfully"}


def install_plugin_dependencies(*, db: Session, current_user: models.User, plugin_id: int) -> dict:
    """Install a plugin's dependencies (admin only; 404 if the plugin is missing)."""
    _require_admin(current_user)

    success = crud_admin.install_plugin_dependencies(db, plugin_id)
    if not success:
        raise NotFoundError("Plugin not found")
    return {"message": "Plugin dependencies installed successfully"}


# =============================================================================
# Plugin Synchronization Management
# =============================================================================

def get_plugin_sync_status(*, current_user: models.User):
    """Return current plugin synchronization status (admin only)."""
    _require_admin(current_user)

    try:
        sync_manager = PluginSyncManager()
        status = sync_manager.get_sync_status()
        return status
    except Exception as e:
        logger.error(f"Failed to get sync status: {e}")
        raise ExternalServiceError(f"Failed to get sync status: {str(e)}", status_code=500)


def sync_plugins_from_repository(*, current_user: models.User):
    """Manually trigger plugin synchronization from repository to database (admin only).

    NOTE (pinned): the inner ``result["success"] is False`` raise is kept as an
    ``HTTPException`` because the bare ``except Exception`` below stringifies the
    caught exception into the final detail (``"Synchronization failed: 500: ..."``);
    converting it would change that byte-identical error body.
    """
    _require_admin(current_user)

    try:
        sync_manager = PluginSyncManager()
        result = sync_manager.sync_plugins_to_database()

        if result["success"]:
            return result
        else:
            raise HTTPException(status_code=500, detail=result.get("error", "Synchronization failed"))

    except Exception as e:
        logger.error(f"Plugin sync failed: {e}")
        raise HTTPException(status_code=500, detail=f"Synchronization failed: {str(e)}")


def get_plugin_consistency_report(*, current_user: models.User, format: str):
    """Return a plugin version consistency report (admin only; json or text)."""
    _require_admin(current_user)

    try:
        validator = PluginVersionValidator()

        if format.lower() == "text":
            # Return plain text report
            report = validator.generate_consistency_report()
            return {"report": report}
        else:
            # Return JSON format
            result = validator.validate_consistency()
            return result

    except Exception as e:
        logger.error(f"Failed to generate consistency report: {e}")
        raise ExternalServiceError(f"Failed to generate report: {str(e)}", status_code=500)


def get_available_plugin_branches(*, current_user: models.User):
    """List available plugin repository branches + the current one (admin only)."""
    _require_admin(current_user)

    try:
        sync_manager = PluginSyncManager()
        branches = sync_manager.get_available_branches()
        current_branch = sync_manager.get_current_branch()

        return {
            "current_branch": current_branch,
            "available_branches": branches,
            "total_branches": len(branches)
        }

    except Exception as e:
        logger.error(f"Failed to get available branches: {e}")
        raise ExternalServiceError(f"Failed to get branches: {str(e)}", status_code=500)


def switch_plugin_branch(*, current_user: models.User, branch: str):
    """Switch the plugin repo branch and re-synchronize the database (admin only).

    NOTE (pinned): inner 400/500 raises are kept as ``HTTPException`` and the
    ``except HTTPException: raise`` guard preserves them; the terminal handler is
    kept verbatim so the swallow-and-stringify error body is byte-identical.
    """
    _require_admin(current_user)

    try:
        sync_manager = PluginSyncManager()

        # Switch branch
        success = sync_manager.switch_branch(branch)
        if not success:
            raise HTTPException(status_code=400, detail=f"Failed to switch to branch: {branch}")

        # Synchronize database
        sync_result = sync_manager.sync_plugins_to_database()

        if sync_result["success"]:
            return {
                "message": f"Successfully switched to branch {branch} and synchronized database",
                "branch": branch,
                "bundle_version": sync_result.get("bundle_version"),
                "sync_result": sync_result
            }
        else:
            raise HTTPException(
                status_code=500,
                detail=f"Branch switched but sync failed: {sync_result.get('error')}"
            )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to switch branch and sync: {e}")
        raise HTTPException(status_code=500, detail=f"Operation failed: {str(e)}")


def quick_consistency_check(*, current_user: models.User):
    """Return a boolean plugin consistency flag with a message (admin only)."""
    _require_admin(current_user)

    try:
        validator = PluginVersionValidator()
        consistent = validator.quick_check()

        return {
            "consistent": consistent,
            "message": "All components are in sync" if consistent else "Inconsistencies detected"
        }

    except Exception as e:
        logger.error(f"Quick consistency check failed: {e}")
        raise ExternalServiceError(f"Check failed: {str(e)}", status_code=500)


# =============================================================================
# Container Manager (메모리 누수 방지)
# =============================================================================

def get_container_manager_status(*, current_user: models.User):
    """Return container-manager tracking status (admin only)."""
    _require_admin(current_user)

    try:
        status = container_manager.get_status()
        return status
    except Exception as e:
        logger.error(f"Failed to get container manager status: {e}")
        raise ExternalServiceError(f"Failed to get status: {str(e)}", status_code=500)


def cleanup_container_manager(*, current_user: models.User):
    """Clean up stale container-manager mappings (admin only).

    NOTE (pinned): the ``"error" in result`` raise is kept as ``HTTPException``
    and guarded by ``except HTTPException: raise`` so its body is byte-identical.
    """
    _require_admin(current_user)

    try:
        result = container_manager.cleanup_stale_mappings()

        if "error" in result:
            raise HTTPException(status_code=500, detail=result["error"])

        logger.info(f"Container manager cleanup completed: {result}")
        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Container manager cleanup failed: {e}")
        raise HTTPException(status_code=500, detail=f"Cleanup failed: {str(e)}")


def force_clear_container_manager(*, current_user: models.User):
    """Force-clear all container-manager mappings (admin only, emergency use)."""
    _require_admin(current_user)

    try:
        result = container_manager.force_clear_all_mappings()
        logger.warning(f"Container manager force cleared: {result}")
        return result
    except Exception as e:
        logger.error(f"Container manager force clear failed: {e}")
        raise ExternalServiceError(f"Force clear failed: {str(e)}", status_code=500)
