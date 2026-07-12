#!/usr/bin/env python
"""
H5AD function-level benchmarks — Before Redis caching baseline.

Standalone script (no pytest) to avoid root conftest DB dependency.

Usage (inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_function_benchmarks.py
"""
import gc
import json
import os
import sys
import tempfile
import time
import tracemalloc
import statistics

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Add backend root to path so we can import app modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from app.datatable.h5ad import (
    get_annotation_columns,
    get_pseudotime_columns,
    organize_column_dtypes,
)

# ---------------------------------------------------------------------------
# Inline benchmark utilities (no relative import needed)
# ---------------------------------------------------------------------------
ITERATIONS = 10


class BenchmarkTimer:
    def __init__(self, name: str, iterations: int = ITERATIONS):
        self.name = name
        self.iterations = iterations
        self.times = []

    def run(self, func, *args, **kwargs):
        result = None
        for _ in range(self.iterations):
            start = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            self.times.append(elapsed)
        return result

    def stats(self):
        if not self.times:
            return {}
        s = sorted(self.times)
        n = len(s)
        return {
            "name": self.name,
            "iterations": n,
            "mean_s": round(statistics.mean(s), 6),
            "median_s": round(statistics.median(s), 6),
            "std_s": round(statistics.stdev(s), 6) if n > 1 else 0.0,
            "min_s": round(s[0], 6),
            "max_s": round(s[-1], 6),
            "p95_s": round(s[int(n * 0.95)], 6) if n >= 20 else round(s[-1], 6),
            "mean_ms": round(statistics.mean(s) * 1000, 2),
            "median_ms": round(statistics.median(s) * 1000, 2),
        }


def measure_memory(func, *args, **kwargs):
    """Run func once, return tracemalloc peak and RSS delta."""
    rss_before = _get_rss_kb()
    tracemalloc.start()
    result = func(*args, **kwargs)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    rss_after = _get_rss_kb()
    return {
        "result": result,
        "tracemalloc_peak_mb": round(peak / (1024 * 1024), 2),
        "rss_before_mb": round(rss_before / 1024, 2),
        "rss_after_mb": round(rss_after / 1024, 2),
        "rss_delta_mb": round((rss_after - rss_before) / 1024, 2),
    }


def _get_rss_kb():
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1])
    except (FileNotFoundError, ValueError):
        pass
    try:
        import resource
        return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except Exception:
        return 0.0


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------
REAL_H5AD = None
for p in ["/app/tutorials/pbmc_light_1000.h5ad",
          os.path.join(os.path.dirname(__file__), "..", "..", "tutorials", "pbmc_light_1000.h5ad")]:
    if os.path.isfile(p):
        REAL_H5AD = os.path.abspath(p)
        break


def create_h5ad(n_obs, n_vars, path):
    adata = ad.AnnData(
        X=np.random.rand(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": np.random.choice(["TypeA", "TypeB", "TypeC", "TypeD"], n_obs),
            "cluster": np.random.choice([f"C{i}" for i in range(10)], n_obs),
            "pseudotime": np.random.rand(n_obs),
        }),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )
    adata.obsm["X_umap"] = np.random.rand(n_obs, 2).astype(np.float32)
    adata.write_h5ad(path)
    size_mb = os.path.getsize(path) / (1024 * 1024)
    del adata
    gc.collect()
    return size_mb


def load_obs(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    del adata
    gc.collect()
    return obs


# ---------------------------------------------------------------------------
# Pipeline simulations (mirror endpoint logic)
# ---------------------------------------------------------------------------
def columns_pipeline(filepath):
    """Mirrors files.py h5ad_columns endpoint."""
    adata = sc.read_h5ad(filepath)
    adata.obs = organize_column_dtypes(adata.obs)
    anno = get_annotation_columns(adata.obs)
    pseudo = get_pseudotime_columns(adata.obs)
    result = {"anno_columns": anno, "pseudo_columns": pseudo}
    del adata
    gc.collect()
    return result


def clusters_pipeline(filepath, anno_column):
    """Mirrors files.py h5ad_cluster endpoint."""
    adata = sc.read_h5ad(filepath)
    adata.obs = organize_column_dtypes(adata.obs)
    clusters = list(map(str, adata.obs[anno_column].value_counts().index))
    result = {"clusters": clusters}
    del adata
    gc.collect()
    return result


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------
def run_bench(name, func, *args, with_memory=True, **kwargs):
    timer = BenchmarkTimer(name, ITERATIONS)
    timer.run(func, *args, **kwargs)
    stats = timer.stats()

    mem = {}
    if with_memory:
        mem_result = measure_memory(func, *args, **kwargs)
        mem = {k: v for k, v in mem_result.items() if k != "result"}

    # Print formatted
    print(f"\n{'='*60}")
    print(f"  {stats['name']}")
    print(f"{'='*60}")
    print(f"  Iterations : {stats['iterations']}")
    print(f"  Mean       : {stats['mean_ms']:.2f} ms")
    print(f"  Median     : {stats['median_ms']:.2f} ms")
    print(f"  Std        : {stats.get('std_s', 0)*1000:.2f} ms")
    print(f"  Min        : {stats['min_s']*1000:.2f} ms")
    print(f"  Max        : {stats['max_s']*1000:.2f} ms")
    if mem:
        print(f"  Peak Mem   : {mem.get('tracemalloc_peak_mb', 'N/A')} MB (tracemalloc)")
        print(f"  RSS Delta  : {mem.get('rss_delta_mb', 'N/A')} MB")

    return {**stats, **mem}


def main():
    all_results = {}

    print("\n" + "#" * 60)
    print("# H5AD Function Benchmark — Before Redis Caching Baseline")
    print("#" * 60)

    # --- Create test data ---
    tmpdir = tempfile.mkdtemp(prefix="bench_h5ad_")
    small_path = os.path.join(tmpdir, "small.h5ad")
    large_path = os.path.join(tmpdir, "large.h5ad")

    print("\n--- Preparing test data ---")
    small_size = create_h5ad(100, 50, small_path)
    print(f"  Small (100×50)   : {small_size:.2f} MB → {small_path}")
    large_size = create_h5ad(10_000, 2_000, large_path)
    print(f"  Large (10k×2k)   : {large_size:.2f} MB → {large_path}")
    if REAL_H5AD:
        real_size = os.path.getsize(REAL_H5AD) / (1024 * 1024)
        print(f"  Real (pbmc 1k)   : {real_size:.2f} MB → {REAL_H5AD}")
    else:
        print("  Real (pbmc 1k)   : NOT FOUND — skipping real data tests")

    datasets = [
        ("small", small_path),
        ("large", large_path),
    ]
    if REAL_H5AD:
        datasets.append(("real", REAL_H5AD))

    # ==================================================================
    # 1. scanpy.read_h5ad
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 1: scanpy.read_h5ad")
    print("=" * 60)

    for label, path in datasets:
        def _read(p=path):
            a = sc.read_h5ad(p)
            shape = a.shape
            del a
            gc.collect()
            return shape

        r = run_bench(f"scanpy.read_h5ad [{label}]", _read)
        all_results[f"read_h5ad_{label}"] = r

    # ==================================================================
    # 2. organize_column_dtypes
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 2: organize_column_dtypes")
    print("=" * 60)

    for label, path in datasets:
        obs = load_obs(path)
        r = run_bench(f"organize_column_dtypes [{label}]", organize_column_dtypes, obs)
        all_results[f"organize_dtypes_{label}"] = r

    # ==================================================================
    # 3. get_annotation_columns
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 3: get_annotation_columns")
    print("=" * 60)

    for label, path in datasets:
        obs = load_obs(path)
        r = run_bench(f"get_annotation_columns [{label}]", get_annotation_columns, obs)
        all_results[f"anno_columns_{label}"] = r

    # ==================================================================
    # 4. get_pseudotime_columns
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 4: get_pseudotime_columns")
    print("=" * 60)

    for label, path in datasets:
        obs = load_obs(path)
        r = run_bench(f"get_pseudotime_columns [{label}]", get_pseudotime_columns, obs)
        all_results[f"pseudo_columns_{label}"] = r

    # ==================================================================
    # 5. Full /columns pipeline
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 5: /columns Full Pipeline (read + organize + anno + pseudo)")
    print("=" * 60)

    for label, path in datasets:
        r = run_bench(f"columns_pipeline [{label}]", columns_pipeline, path)
        all_results[f"columns_pipeline_{label}"] = r

    # ==================================================================
    # 6. Full /clusters pipeline
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 6: /clusters Full Pipeline (read + organize + value_counts)")
    print("=" * 60)

    for label, path in datasets:
        # Determine a valid annotation column
        if label == "real":
            obs = load_obs(path)
            anno_cols = get_annotation_columns(obs)
            anno_col = anno_cols[0] if anno_cols else "cluster"
        else:
            anno_col = "cluster"

        r = run_bench(f"clusters_pipeline [{label}]", clusters_pipeline, path, anno_col)
        all_results[f"clusters_pipeline_{label}"] = r

    # ==================================================================
    # Summary table
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SUMMARY TABLE")
    print("=" * 60)
    print(f"{'Benchmark':<45} {'Mean (ms)':>10} {'Median (ms)':>12} {'Peak Mem (MB)':>14}")
    print("-" * 85)
    for key, val in all_results.items():
        mem_str = str(val.get("tracemalloc_peak_mb", "N/A"))
        print(f"{val['name']:<45} {val['mean_ms']:>10.2f} {val['median_ms']:>12.2f} {mem_str:>14}")

    # ==================================================================
    # JSON output
    # ==================================================================
    json_path = os.path.join(tmpdir, "benchmark_results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Raw JSON saved → {json_path}")

    # Also print JSON to stdout for capture
    print("\n--- RAW JSON ---")
    print(json.dumps(all_results, indent=2))

    # Cleanup
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

    return all_results


if __name__ == "__main__":
    main()
