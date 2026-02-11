"""
DataTable DataFrame caching via Redis.

Caches the parsed DataFrame (obs + UMAP/PCA) from load_tab_file()
so that page/sort/filter changes don't re-read the entire h5ad file.

Cache key includes file mtime for automatic invalidation on file changes.
Pickle is used for serialization — safe because only trusted internal
DataFrames are serialized, never external/user-supplied data.

SECURITY: Never store user-supplied or externally-sourced objects in this cache.
pickle.loads() can execute arbitrary code — only cache internally-generated DataFrames.
"""
import os
import pickle
import logging
from typing import Optional

import pandas as pd

from app.common.utils.redis_cache import cache_get_bytes, cache_set_bytes, cache_delete_pattern
from app.common.config import settings

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
