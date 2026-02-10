"""
h5ad metadata in-process cache using cachetools.TTLCache.

Caches lightweight metadata (column names, cluster lists) to avoid
re-reading entire h5ad files via scanpy.read_h5ad() on every request.

Design decisions:
- In-process TTLCache (not Redis) because cached data is <1KB per entry,
  backend runs a single Uvicorn worker, and Celery never calls these endpoints.
- Cache key includes file mtime for automatic invalidation on file changes.
- Thread lock guards concurrent access from async endpoint threads.
- Explicit invalidate_file() called on file deletion for immediate purge.
"""

import os
import threading
from typing import Optional, Dict, Any, List

from cachetools import TTLCache

_lock = threading.Lock()
_columns_cache: TTLCache = TTLCache(maxsize=500, ttl=86400)  # 24h
_clusters_cache: TTLCache = TTLCache(maxsize=500, ttl=86400)


def _make_columns_key(filepath: str) -> str:
    mtime = os.path.getmtime(filepath)
    return f"{filepath}:{mtime}"


def _make_clusters_key(filepath: str, anno_column: str) -> str:
    mtime = os.path.getmtime(filepath)
    return f"{filepath}:{anno_column}:{mtime}"


def get_cached_columns(filepath: str) -> Optional[Dict[str, Any]]:
    key = _make_columns_key(filepath)
    with _lock:
        return _columns_cache.get(key)


def set_cached_columns(filepath: str, data: Dict[str, Any]) -> None:
    key = _make_columns_key(filepath)
    with _lock:
        _columns_cache[key] = data


def get_cached_clusters(filepath: str, anno_column: str) -> Optional[List[str]]:
    key = _make_clusters_key(filepath, anno_column)
    with _lock:
        return _clusters_cache.get(key)


def set_cached_clusters(filepath: str, anno_column: str, data: List[str]) -> None:
    key = _make_clusters_key(filepath, anno_column)
    with _lock:
        _clusters_cache[key] = data


def invalidate_file(filepath: str) -> int:
    """Remove all cache entries for the given file path."""
    count = 0
    with _lock:
        for cache in (_columns_cache, _clusters_cache):
            keys_to_del = [k for k in cache if k.startswith(filepath + ":")]
            for k in keys_to_del:
                del cache[k]
                count += 1
    return count
