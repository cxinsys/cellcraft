"""
Plugin list caching via Redis.

Caches the /plugin/list response per user so that repeated modal opens
(Algorithm, Visualization) don't re-query the database and Docker daemon.

Uses JSON serialization — plugin list data is plain dicts/strings.
On Redis failure, operations silently return None / skip,
allowing the application to fall back to the original query path.
"""
import json
import logging
from typing import Optional

from app.common.utils.redis_cache import cache_get_bytes, cache_set_bytes, cache_delete_pattern
from app.common.config import settings

logger = logging.getLogger(__name__)


def _make_key(user_id: int) -> str:
    return f"plugin:list:{user_id}"


def get_cached_plugin_list(user_id: int) -> Optional[dict]:
    """Get cached plugin list response for a user. Returns None on miss or error."""
    data = cache_get_bytes(_make_key(user_id))
    if data is None:
        return None
    try:
        return json.loads(data)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        logger.warning(f"Plugin cache JSON decode failed for user {user_id}: {e}")
        return None


def set_cached_plugin_list(user_id: int, data: dict) -> bool:
    """Cache plugin list response for a user."""
    try:
        payload = json.dumps(data, default=str).encode("utf-8")
        return cache_set_bytes(_make_key(user_id), payload, settings.PLUGIN_LIST_CACHE_TTL)
    except (TypeError, ValueError) as e:
        logger.warning(f"Plugin cache JSON encode failed for user {user_id}: {e}")
        return False


def invalidate_plugin_cache_for_user(user_id: int) -> int:
    """Invalidate cached plugin list for a specific user."""
    return cache_delete_pattern(f"plugin:list:{user_id}")


def invalidate_all_plugin_cache() -> int:
    """Invalidate cached plugin list for all users."""
    return cache_delete_pattern("plugin:list:*")
