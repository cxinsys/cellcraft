from celery.result import AsyncResult
from celery.signals import worker_ready
import logging

# create_celery moved to app.worker.app in PR-2 (phase-1-app-assembly).
# Re-exported here as a shim to preserve the existing import path; the shim is
# removed in Phase 2 once callers migrate to app.worker.app.
from app.worker.app import create_celery  # noqa: F401

logger = logging.getLogger(__name__)


@worker_ready.connect
def on_worker_ready(sender, **kwargs):
    """Initialize resource totals and clean up stale allocations on worker startup."""
    from app.shared.resources import (
        initialize_resource_totals, cleanup_stale_resources
    )
    result = initialize_resource_totals()
    cleanup_stale_resources()
    logger.info(
        f"Worker ready — resources initialized: "
        f"CPU={result['cpu_total']}, GPU={result['gpu_total']}, "
        f"concurrency={result['worker_concurrency']}"
    )


def get_task_info(task_id):
    """
    return task info for the given task_id
    """
    task_result = AsyncResult(task_id)
    result = {
        "task_id": task_id,
        "task_status": task_result.status,
        "task_result": task_result.result,
    }
    return result
