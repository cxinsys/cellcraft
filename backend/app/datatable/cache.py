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


# ---------------------------------------------------------------------------
# DataTable DataFrame caching via Redis.
#
# Caches the parsed DataFrame (obs + UMAP/PCA) from load_tab_file()
# so that page/sort/filter changes don't re-read the entire h5ad file.
#
# Cache key includes file mtime for automatic invalidation on file changes.
# Pickle is used for serialization — safe because only trusted internal
# DataFrames are serialized, never external/user-supplied data.
#
# SECURITY: Never store user-supplied or externally-sourced objects in this cache.
# pickle.loads() can execute arbitrary code — only cache internally-generated DataFrames.
#
# Merged from app/common/utils/datatable_cache.py in PR-4 (phase-2-domain-skeleton).
# ---------------------------------------------------------------------------
import pickle
import logging

import pandas as pd

from app.shared.redis import cache_get_bytes, cache_set_bytes, cache_delete_pattern
from app.core.config import settings

logger = logging.getLogger(__name__)

MAX_CACHEABLE_SIZE_MB = 200


def _make_key(filepath: str) -> str:
    mtime = os.path.getmtime(filepath)
    return f"datatable:df:{filepath}:{mtime}"


def get_cached_dataframe(filepath: str) -> Optional[pd.DataFrame]:
    key = _make_key(filepath)
    data = cache_get_bytes(key)
    if data is None:
        return None
    try:
        return pickle.loads(data)
    except Exception as e:
        logger.warning(f"DataFrame unpickle failed for {key}: {e}")
        return None


def set_cached_dataframe(filepath: str, df: pd.DataFrame) -> bool:
    size_mb = df.memory_usage(deep=True).sum() / (1024 * 1024)
    if size_mb > MAX_CACHEABLE_SIZE_MB:
        logger.info(f"DataFrame too large to cache: {size_mb:.0f}MB > {MAX_CACHEABLE_SIZE_MB}MB")
        return False

    key = _make_key(filepath)
    try:
        data = pickle.dumps(df, protocol=pickle.HIGHEST_PROTOCOL)
        return cache_set_bytes(key, data, settings.DATATABLE_DF_CACHE_TTL)
    except Exception as e:
        logger.warning(f"DataFrame pickle/cache failed for {key}: {e}")
        return False


def invalidate_datatable(filepath: str) -> int:
    return cache_delete_pattern(f"datatable:df:{filepath}:*")
