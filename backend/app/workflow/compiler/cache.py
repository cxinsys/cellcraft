"""DAG/로그 파일 상태 캐시 (PR-10 split from ``dag_parser``).

``SnakemakeDAGParser`` (parsing) 및 ``SnakemakeRuleStatusTracker`` (status)에서
공유하는 캐시 인스턴스(``_dag_cache``, ``_log_cache``)와 캐시 관리 유틸리티를
담는다. 기존 ``dag_parser`` 모듈에 정의돼 있던 것과 동일하며, 파사드가 그대로
re-export 한다.
"""
import os
import time
import hashlib
import logging
import threading
from typing import Dict, Any, Optional, Tuple
from functools import lru_cache

# 로거 설정
logger = logging.getLogger(__name__)


class DAGCache:
    """DAG 구조와 로그 상태를 캐싱하는 클래스"""
    def __init__(self, max_size: int = 100, ttl_seconds: int = 300):
        self.max_size = max_size
        self.ttl_seconds = ttl_seconds  # 5분 TTL
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._access_times: Dict[str, float] = {}
        self._file_mtimes: Dict[str, float] = {}
        self._lock = threading.Lock()

    def _get_cache_key(self, file_path: str) -> str:
        """파일 경로와 수정 시간을 기반으로 캐시 키 생성"""
        try:
            file_stat = os.stat(file_path)
            mtime = file_stat.st_mtime
            size = file_stat.st_size

            # 파일 경로, 수정시간, 크기를 해시화
            key_data = f"{file_path}_{mtime}_{size}"
            return hashlib.md5(key_data.encode()).hexdigest()
        except (OSError, FileNotFoundError):
            # 파일이 없거나 접근할 수 없는 경우
            return hashlib.md5(file_path.encode()).hexdigest()

    def _is_cache_valid(self, key: str, file_path: str) -> bool:
        """캐시 유효성 검사"""
        try:
            current_time = time.time()

            # TTL 확인
            if key not in self._access_times:
                return False

            if current_time - self._access_times[key] > self.ttl_seconds:
                return False

            # 파일 수정 시간 확인
            if os.path.exists(file_path):
                current_mtime = os.path.getmtime(file_path)
                cached_mtime = self._file_mtimes.get(key)

                if cached_mtime is None or current_mtime > cached_mtime:
                    return False

            return True

        except Exception as e:
            logger.warning(f"Error checking cache validity for {file_path}: {e}")
            return False

    def get(self, file_path: str) -> Optional[Dict[str, Any]]:
        """캐시에서 DAG 데이터 조회"""
        with self._lock:
            key = self._get_cache_key(file_path)

            if key in self._cache and self._is_cache_valid(key, file_path):
                # 액세스 시간 업데이트
                self._access_times[key] = time.time()
                logger.debug(f"Cache hit for {file_path}")
                return self._cache[key].copy()

            # 캐시 미스 또는 무효한 캐시
            if key in self._cache:
                logger.debug(f"Cache invalidated for {file_path}")
                self._remove_key(key)

            return None

    def set(self, file_path: str, data: Dict[str, Any]) -> None:
        """DAG 데이터를 캐시에 저장"""
        with self._lock:
            key = self._get_cache_key(file_path)
            current_time = time.time()

            # 캐시 크기 제한 확인
            if len(self._cache) >= self.max_size:
                self._evict_oldest()

            # 파일 수정 시간 저장
            try:
                if os.path.exists(file_path):
                    self._file_mtimes[key] = os.path.getmtime(file_path)
            except Exception as e:
                logger.warning(f"Error getting mtime for {file_path}: {e}")

            self._cache[key] = data.copy()
            self._access_times[key] = current_time
            logger.debug(f"Cached DAG data for {file_path}")

    def _remove_key(self, key: str) -> None:
        """캐시에서 키 제거"""
        self._cache.pop(key, None)
        self._access_times.pop(key, None)
        self._file_mtimes.pop(key, None)

    def _evict_oldest(self) -> None:
        """가장 오래된 캐시 엔트리 제거"""
        if not self._access_times:
            return

        oldest_key = min(self._access_times.keys(), key=lambda k: self._access_times[k])
        self._remove_key(oldest_key)
        logger.debug(f"Evicted oldest cache entry: {oldest_key}")

    def _evict_batch(self, count: int = 10) -> int:
        """
        여러 캐시 엔트리를 한 번에 제거 (메모리 압박 시 효율적 정리)

        Args:
            count: 제거할 엔트리 수 (기본값: 10)

        Returns:
            int: 실제로 제거된 엔트리 수
        """
        evicted = 0
        for _ in range(min(count, len(self._access_times))):
            if self._access_times:
                self._evict_oldest()
                evicted += 1
        if evicted > 0:
            logger.debug(f"Batch evicted {evicted} cache entries")
        return evicted

    def clear(self) -> None:
        """전체 캐시 초기화"""
        with self._lock:
            self._cache.clear()
            self._access_times.clear()
            self._file_mtimes.clear()
            logger.info("DAG cache cleared")

    def clear_all(self) -> int:
        """
        전체 캐시 정리 (관리자용) - 제거된 엔트리 수 반환

        Returns:
            int: 정리된 캐시 엔트리 수
        """
        with self._lock:
            count = len(self._cache)
            self._cache.clear()
            self._access_times.clear()
            self._file_mtimes.clear()
            logger.info(f"DAG cache cleared: {count} entries removed")
            return count

    def get_stats(self) -> Dict[str, Any]:
        """캐시 통계 반환"""
        with self._lock:
            return {
                "size": len(self._cache),
                "max_size": self.max_size,
                "ttl_seconds": self.ttl_seconds,
                "file_mtimes_count": len(self._file_mtimes),
                "access_times_count": len(self._access_times)
            }

# 전역 캐시 인스턴스
_dag_cache = DAGCache()


class LogFileCache:
    """로그 파일 상태를 캐싱하는 클래스"""
    def __init__(self, ttl_seconds: int = 60):  # 1분 TTL
        self.ttl_seconds = ttl_seconds
        self._cache: Dict[str, Tuple[float, bool]] = {}  # path -> (timestamp, exists)
        self._lock = threading.Lock()

    def file_exists(self, file_path: str) -> bool:
        """파일 존재 여부를 캐시를 통해 확인"""
        with self._lock:
            current_time = time.time()

            if file_path in self._cache:
                cached_time, cached_exists = self._cache[file_path]

                # TTL 확인
                if current_time - cached_time < self.ttl_seconds:
                    return cached_exists

            # 캐시 미스 또는 만료
            try:
                exists = os.path.exists(file_path)
                self._cache[file_path] = (current_time, exists)
                return exists
            except Exception as e:
                logger.warning(f"Error checking file existence {file_path}: {e}")
                return False

    def clear(self) -> None:
        """캐시 초기화"""
        with self._lock:
            self._cache.clear()

# 전역 로그 캐시 인스턴스
_log_cache = LogFileCache()


# 캐시 관리 유틸리티 함수들
def clear_all_caches():
    """모든 캐시를 초기화"""
    _dag_cache.clear()
    _log_cache.clear()
    logger.info("All caches cleared")

def get_cache_stats() -> Dict[str, Any]:
    """캐시 통계 정보 반환"""
    return {
        "dag_cache": _dag_cache.get_stats(),
        "log_cache": {
            "size": len(_log_cache._cache),
            "ttl_seconds": _log_cache.ttl_seconds
        }
    }

def set_cache_config(dag_max_size: int = None, dag_ttl: int = None, log_ttl: int = None):
    """캐시 설정 변경"""
    global _dag_cache, _log_cache

    if dag_max_size is not None:
        _dag_cache.max_size = dag_max_size

    if dag_ttl is not None:
        _dag_cache.ttl_seconds = dag_ttl

    if log_ttl is not None:
        _log_cache.ttl_seconds = log_ttl

    logger.info(f"Cache config updated: DAG(max_size={_dag_cache.max_size}, ttl={_dag_cache.ttl_seconds}s), Log(ttl={_log_cache.ttl_seconds}s)")

@lru_cache(maxsize=128)
def get_file_hash(file_path: str, mtime: float) -> str:
    """파일 해시 계산 (LRU 캐시 적용)"""
    try:
        hash_data = f"{file_path}_{mtime}"
        return hashlib.md5(hash_data.encode()).hexdigest()
    except Exception as e:
        logger.error(f"Error calculating file hash for {file_path}: {e}")
        return hashlib.md5(file_path.encode()).hexdigest()
