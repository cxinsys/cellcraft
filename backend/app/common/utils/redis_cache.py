"""
Redis cache layer with graceful fallback.

On Redis failure, operations silently return None (get) or skip (set),
allowing the application to fall back to the original data loading path.
"""
import logging
from typing import Optional

import redis

from app.common.config import settings

logger = logging.getLogger(__name__)

_pool: Optional[redis.ConnectionPool] = None


def _get_pool() -> redis.ConnectionPool:
    global _pool
    if _pool is None:
        _pool = redis.ConnectionPool.from_url(
            settings.REDIS_URL,
            max_connections=20,
            decode_responses=False,
        )
    return _pool


def get_redis_client() -> redis.Redis:
    return redis.Redis(connection_pool=_get_pool())


def cache_get_bytes(key: str) -> Optional[bytes]:
    """Get binary data from cache (for DataFrame serialization)."""
    try:
        client = get_redis_client()
        return client.get(key)
    except redis.RedisError as e:
        logger.warning(f"Redis cache_get_bytes failed for {key}: {e}")
        return None


def cache_set_bytes(key: str, value: bytes, ttl: int) -> bool:
    """Set binary data in cache with TTL."""
    try:
        client = get_redis_client()
        client.setex(key, ttl, value)
        return True
    except redis.RedisError as e:
        logger.warning(f"Redis cache_set_bytes failed for {key}: {e}")
        return False


def cache_delete_pattern(pattern: str) -> int:
    """Delete all keys matching a pattern using SCAN (non-blocking)."""
    try:
        client = get_redis_client()
        count = 0
        for key in client.scan_iter(match=pattern, count=100):
            client.delete(key)
            count += 1
        return count
    except redis.RedisError as e:
        logger.warning(f"Redis cache_delete_pattern failed for {pattern}: {e}")
        return 0
