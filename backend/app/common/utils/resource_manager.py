"""
Redis-based resource semaphore for Celery task queue gating.

Uses Lua scripts for atomic check-and-increment operations to prevent
race conditions. On Redis failure, all operations fall back to ungated
execution (acquire returns True, release is a no-op).

Redis key schema:
  resource:cpu:total   - Total CPU slots
  resource:gpu:total   - Total GPU slots
  resource:cpu:used    - Currently used CPU slots
  resource:gpu:used    - Currently used GPU slots
  resource:task:{id}   - Per-task allocation hash (type, slots, plugin_name, started_at)
"""
import logging
import os
from datetime import datetime, timezone
from typing import Optional

import redis

from app.common.config import settings
from app.common.utils.redis_cache import get_redis_client

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lua scripts for atomic Redis operations
# ---------------------------------------------------------------------------

# Atomically check available capacity and increment used counter.
# Returns 1 on success (slots acquired), 0 on failure (insufficient capacity).
_LUA_ACQUIRE = """
local total_key = KEYS[1]
local used_key  = KEYS[2]
local task_key  = KEYS[3]
local requested = tonumber(ARGV[1])
local rtype     = ARGV[2]
local plugin    = ARGV[3]
local started   = ARGV[4]

local total = tonumber(redis.call('GET', total_key) or 0)
local used  = tonumber(redis.call('GET', used_key) or 0)

if used + requested <= total then
    redis.call('INCRBY', used_key, requested)
    redis.call('HSET', task_key,
        'type', rtype,
        'slots', requested,
        'plugin_name', plugin,
        'started_at', started)
    return 1
end
return 0
"""

# Atomically decrement used counter and delete task key.
# Ensures used never goes below 0.
_LUA_RELEASE = """
local used_key = KEYS[1]
local task_key = KEYS[2]
local slots    = tonumber(ARGV[1])

local used = tonumber(redis.call('GET', used_key) or 0)
local new_used = used - slots
if new_used < 0 then
    new_used = 0
end
redis.call('SET', used_key, new_used)
redis.call('DEL', task_key)
return new_used
"""


# ---------------------------------------------------------------------------
# Auto-detection helpers
# ---------------------------------------------------------------------------

def detect_cpu_count() -> int:
    """Detect total CPU slots. Uses config override if set, else os.cpu_count()."""
    if settings.RESOURCE_CPU_TOTAL > 0:
        return settings.RESOURCE_CPU_TOTAL
    count = os.cpu_count()
    return count if count else 4  # fallback


def detect_gpu_count() -> int:
    """Detect total GPU devices. Priority: config > GPU_COUNT env > GPUtil > 0."""
    if settings.RESOURCE_GPU_TOTAL > 0:
        return settings.RESOURCE_GPU_TOTAL
    env_val = os.environ.get('GPU_COUNT')
    if env_val is not None:
        try:
            return int(env_val)
        except (ValueError, TypeError):
            pass
    try:
        import GPUtil
        return len(GPUtil.getGPUs())
    except Exception:
        pass
    return 0


def detect_worker_concurrency(cpu_total: int) -> int:
    """Calculate optimal worker concurrency from CPU total and default slot size."""
    if settings.CELERY_WORKER_CONCURRENCY > 0:
        return settings.CELERY_WORKER_CONCURRENCY
    return max(1, cpu_total // settings.RESOURCE_DEFAULT_CPU_SLOTS)


# ---------------------------------------------------------------------------
# Core resource operations
# ---------------------------------------------------------------------------

def acquire_slots(resource_type: str, slots: int, task_id: str,
                  plugin_name: str = '') -> bool:
    """
    Atomically acquire resource slots for a task.

    Returns True if slots were acquired, False if insufficient capacity.
    On Redis failure, returns True (graceful fallback — ungated execution).
    """
    try:
        client = get_redis_client()
        total_key = f"resource:{resource_type}:total"
        used_key = f"resource:{resource_type}:used"
        task_key = f"resource:task:{task_id}"
        started = datetime.now(timezone.utc).isoformat()

        result = client.eval(
            _LUA_ACQUIRE, 3,
            total_key, used_key, task_key,
            slots, resource_type, plugin_name, started,
        )
        acquired = int(result) == 1
        if acquired:
            logger.info(
                f"Resource acquired: task={task_id} type={resource_type} "
                f"slots={slots} plugin={plugin_name}"
            )
        else:
            logger.info(
                f"Resource unavailable: task={task_id} type={resource_type} "
                f"slots={slots} — will retry"
            )
        return acquired

    except redis.RedisError as e:
        logger.warning(
            f"Redis error in acquire_slots for task {task_id}: {e}. "
            "Falling back to ungated execution."
        )
        return True


def release_slots(resource_type: str, slots: int, task_id: str) -> None:
    """
    Atomically release resource slots for a task.
    Safe to call multiple times — Lua script prevents negative counters.
    """
    try:
        client = get_redis_client()
        used_key = f"resource:{resource_type}:used"
        task_key = f"resource:task:{task_id}"
        client.eval(_LUA_RELEASE, 2, used_key, task_key, slots)
        logger.info(
            f"Resource released: task={task_id} type={resource_type} slots={slots}"
        )
    except redis.RedisError as e:
        logger.warning(f"Redis error in release_slots for task {task_id}: {e}")


def release_slots_by_task_id(task_id: str) -> bool:
    """
    Look up the task's allocation hash and release its slots.
    Used as a safety net in after_return when kwargs may not be available.
    Returns True if slots were released, False if no allocation found.
    """
    try:
        client = get_redis_client()
        task_key = f"resource:task:{task_id}"
        info = client.hgetall(task_key)
        if not info:
            return False

        resource_type = info.get(b'type', info.get('type', b'cpu'))
        slots = info.get(b'slots', info.get('slots', b'0'))
        # Decode bytes if needed
        if isinstance(resource_type, bytes):
            resource_type = resource_type.decode()
        if isinstance(slots, bytes):
            slots = slots.decode()

        release_slots(resource_type, int(slots), task_id)
        return True

    except redis.RedisError as e:
        logger.warning(f"Redis error in release_slots_by_task_id for {task_id}: {e}")
        return False


# ---------------------------------------------------------------------------
# Initialization (called on worker_ready)
# ---------------------------------------------------------------------------

def initialize_resource_totals() -> dict:
    """
    Detect system resources and write totals to Redis.
    Called once on worker startup via worker_ready signal.
    Returns dict with detected values for logging.
    """
    cpu_total = detect_cpu_count()
    gpu_total = detect_gpu_count()
    concurrency = detect_worker_concurrency(cpu_total)

    cpu_source = "config" if settings.RESOURCE_CPU_TOTAL > 0 else "auto"
    gpu_source = (
        "config" if settings.RESOURCE_GPU_TOTAL > 0
        else ("env:GPU_COUNT" if os.environ.get('GPU_COUNT') else "auto")
    )

    try:
        client = get_redis_client()
        client.set("resource:cpu:total", cpu_total)
        client.set("resource:gpu:total", gpu_total)
        logger.info(
            f"Resource detection: CPU={cpu_total} ({cpu_source}), "
            f"GPU={gpu_total} ({gpu_source}), "
            f"Worker concurrency={concurrency}"
        )
    except redis.RedisError as e:
        logger.warning(f"Redis error in initialize_resource_totals: {e}")

    return {
        'cpu_total': cpu_total,
        'gpu_total': gpu_total,
        'worker_concurrency': concurrency,
    }


def cleanup_stale_resources() -> None:
    """
    Reset used counters and remove stale task keys on worker startup.
    Safe for single-worker architecture — assumes no other worker holds slots.
    """
    try:
        client = get_redis_client()
        # Reset used counters
        client.set("resource:cpu:used", 0)
        client.set("resource:gpu:used", 0)
        # Remove stale task allocation keys
        count = 0
        for key in client.scan_iter(match="resource:task:*", count=100):
            client.delete(key)
            count += 1
        if count > 0:
            logger.info(f"Cleaned up {count} stale resource task key(s)")
    except redis.RedisError as e:
        logger.warning(f"Redis error in cleanup_stale_resources: {e}")


# ---------------------------------------------------------------------------
# Monitoring API helper
# ---------------------------------------------------------------------------

def get_resource_status() -> Optional[dict]:
    """
    Return comprehensive resource status for the monitoring API.
    Returns None on Redis failure (caller should return HTTP 503).
    """
    try:
        client = get_redis_client()
        cpu_total = int(client.get("resource:cpu:total") or 0)
        gpu_total = int(client.get("resource:gpu:total") or 0)
        cpu_used = int(client.get("resource:cpu:used") or 0)
        gpu_used = int(client.get("resource:gpu:used") or 0)

        # Collect running task allocations
        tasks = []
        for key in client.scan_iter(match="resource:task:*", count=100):
            info = client.hgetall(key)
            if not info:
                continue
            # Extract task_id from key pattern "resource:task:{task_id}"
            key_str = key.decode() if isinstance(key, bytes) else key
            task_id = key_str.replace("resource:task:", "", 1)

            def _decode(val):
                return val.decode() if isinstance(val, bytes) else val

            tasks.append({
                'task_id': task_id,
                'resource_type': _decode(info.get(b'type', info.get('type', b''))),
                'resource_slots': int(_decode(info.get(b'slots', info.get('slots', b'0')))),
                'plugin_name': _decode(info.get(b'plugin_name', info.get('plugin_name', b''))),
                'started_at': _decode(info.get(b'started_at', info.get('started_at', b''))),
            })

        return {
            'cpu': {
                'total': cpu_total,
                'used': cpu_used,
                'available': max(0, cpu_total - cpu_used),
            },
            'gpu': {
                'total': gpu_total,
                'used': gpu_used,
                'available': max(0, gpu_total - gpu_used),
            },
            'tasks': tasks,
        }

    except redis.RedisError as e:
        logger.warning(f"Redis error in get_resource_status: {e}")
        return None
