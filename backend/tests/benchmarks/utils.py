"""
Benchmark utility module for h5ad performance measurement.

Provides BenchmarkTimer and MemoryProfiler for systematic
before/after performance comparison of Redis caching.
"""
import time
import tracemalloc
import statistics
import json
import os
from typing import Callable, Any


class BenchmarkTimer:
    """Measures execution time using time.perf_counter() with repeated runs."""

    def __init__(self, name: str, iterations: int = 10):
        self.name = name
        self.iterations = iterations
        self.times: list[float] = []

    def run(self, func: Callable, *args, **kwargs) -> Any:
        """Execute func for self.iterations times, recording each duration."""
        result = None
        for i in range(self.iterations):
            start = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            self.times.append(elapsed)
        return result

    def stats(self) -> dict:
        if not self.times:
            return {}
        sorted_times = sorted(self.times)
        n = len(sorted_times)
        return {
            "name": self.name,
            "iterations": n,
            "mean_s": statistics.mean(sorted_times),
            "median_s": statistics.median(sorted_times),
            "std_s": statistics.stdev(sorted_times) if n > 1 else 0.0,
            "min_s": sorted_times[0],
            "max_s": sorted_times[-1],
            "p95_s": sorted_times[int(n * 0.95)] if n >= 20 else sorted_times[-1],
            "p99_s": sorted_times[int(n * 0.99)] if n >= 100 else sorted_times[-1],
        }


class MemoryProfiler:
    """Measures memory usage via tracemalloc (peak) and /proc RSS."""

    @staticmethod
    def measure(func: Callable, *args, **kwargs) -> dict:
        """Run func once, return peak tracemalloc and RSS delta."""
        # RSS before (Linux /proc/self/status)
        rss_before = MemoryProfiler._get_rss_kb()

        tracemalloc.start()
        result = func(*args, **kwargs)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        rss_after = MemoryProfiler._get_rss_kb()

        return {
            "result": result,
            "tracemalloc_peak_mb": round(peak / (1024 * 1024), 2),
            "rss_before_mb": round(rss_before / 1024, 2),
            "rss_after_mb": round(rss_after / 1024, 2),
            "rss_delta_mb": round((rss_after - rss_before) / 1024, 2),
        }

    @staticmethod
    def _get_rss_kb() -> float:
        """Read VmRSS from /proc/self/status (Linux only)."""
        try:
            with open("/proc/self/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        return float(line.split()[1])  # kB
        except (FileNotFoundError, ValueError):
            pass
        # Fallback: resource module
        try:
            import resource
            return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        except Exception:
            return 0.0


def format_benchmark_result(timer_stats: dict, memory_info: dict | None = None) -> str:
    """Pretty-print a single benchmark result."""
    lines = [
        f"=== {timer_stats['name']} ===",
        f"  Iterations : {timer_stats['iterations']}",
        f"  Mean       : {timer_stats['mean_s']*1000:.2f} ms",
        f"  Median     : {timer_stats['median_s']*1000:.2f} ms",
        f"  Std        : {timer_stats['std_s']*1000:.2f} ms",
        f"  Min        : {timer_stats['min_s']*1000:.2f} ms",
        f"  Max        : {timer_stats['max_s']*1000:.2f} ms",
    ]
    if memory_info:
        lines.extend([
            f"  Peak Mem   : {memory_info.get('tracemalloc_peak_mb', 'N/A')} MB (tracemalloc)",
            f"  RSS Delta  : {memory_info.get('rss_delta_mb', 'N/A')} MB",
        ])
    return "\n".join(lines)


def get_real_h5ad_path() -> str | None:
    """Return path to pbmc_light_1000.h5ad if it exists."""
    candidates = [
        "/app/tutorials/pbmc_light_1000.h5ad",  # inside container
        os.path.join(os.path.dirname(__file__), "..", "..", "tutorials", "pbmc_light_1000.h5ad"),
    ]
    for p in candidates:
        resolved = os.path.abspath(p)
        if os.path.isfile(resolved):
            return resolved
    return None
