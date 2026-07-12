# v1.0.7 리팩토링 shim — PR-4에서 제거
from app.shared.resources import *  # noqa: F401,F403
from app.shared.resources import (  # noqa: F401
    detect_cpu_count,
    detect_gpu_count,
    detect_worker_concurrency,
    acquire_slots,
    release_slots,
    release_slots_by_task_id,
    initialize_resource_totals,
    cleanup_stale_resources,
    get_resource_status,
)
