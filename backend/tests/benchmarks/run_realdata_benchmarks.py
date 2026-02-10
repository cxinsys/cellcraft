#!/usr/bin/env python
"""
Real scRNA-seq data benchmarks — Before Redis caching baseline.

Measures function-level and API-level performance using actual h5ad files
(pbmc_light_1000, CRCBraw, GTEx_8_tissues) to establish realistic baselines.

Usage (inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_realdata_benchmarks.py
"""
import gc
import json
import os
import sys
import time
import statistics
import tracemalloc

import scanpy as sc

# Add backend root to path so we can import app modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from app.common.utils.h5ad_utils import (
    get_annotation_columns,
    get_pseudotime_columns,
    organize_column_dtypes,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATASETS = {
    "pbmc": {
        "path": "/app/tutorials/pbmc_light_1000.h5ad",
        "iterations_func": 10,
        "iterations_api": 10,
        "description": "1K cells × 13.7K genes, 27 MB",
    },
    "CRC": {
        "path": "/app/tutorials/CRCBraw.h5ad",
        "iterations_func": 10,
        "iterations_api": 10,
        "description": "12.6K cells × 27.9K genes, 138 MB",
    },
    "GTEx": {
        "path": "/app/tutorials/GTEx_8_tissues.h5ad",
        "iterations_func": 5,
        "iterations_api": 3,
        "description": "209K cells × 17.7K genes, 1.8 GB",
    },
}

# ---------------------------------------------------------------------------
# Benchmark utilities
# ---------------------------------------------------------------------------

class BenchmarkTimer:
    def __init__(self, name: str, iterations: int):
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
# Helper functions
# ---------------------------------------------------------------------------

def load_obs(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    del adata
    gc.collect()
    return obs


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

def run_bench(name, func, iterations, *args, with_memory=True, **kwargs):
    timer = BenchmarkTimer(name, iterations)
    timer.run(func, *args, **kwargs)
    stats = timer.stats()

    mem = {}
    if with_memory:
        mem_result = measure_memory(func, *args, **kwargs)
        mem = {k: v for k, v in mem_result.items() if k != "result"}

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


def get_data_info(path):
    """Read h5ad once to get shape, obs columns, and sparsity info."""
    adata = sc.read_h5ad(path)
    info = {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "n_obs_columns": len(adata.obs.columns),
        "file_size_mb": round(os.path.getsize(path) / (1024 * 1024), 1),
    }
    # Sparsity
    try:
        from scipy import sparse
        if sparse.issparse(adata.X):
            info["sparse"] = True
            info["dtype"] = str(adata.X.dtype)
            total = adata.X.shape[0] * adata.X.shape[1]
            nnz = adata.X.nnz
            info["sparsity_pct"] = round((1 - nnz / total) * 100, 1)
        else:
            info["sparse"] = False
            info["dtype"] = str(adata.X.dtype)
    except Exception:
        pass

    # Find annotation column
    obs_organized = organize_column_dtypes(adata.obs)
    anno_cols = get_annotation_columns(obs_organized)
    info["anno_columns_detected"] = anno_cols
    info["first_anno_col"] = anno_cols[0] if anno_cols else None

    del adata
    gc.collect()
    return info


# ---------------------------------------------------------------------------
# PART 1: Function-level benchmarks
# ---------------------------------------------------------------------------

def run_function_benchmarks():
    all_results = {}

    print("\n" + "#" * 60)
    print("# PART 1: Function-Level Benchmarks (Real scRNA-seq Data)")
    print("#" * 60)

    # Collect dataset info
    data_info = {}
    for label, cfg in DATASETS.items():
        path = cfg["path"]
        if not os.path.isfile(path):
            print(f"\n  WARNING: {path} not found — skipping {label}")
            continue
        print(f"\n--- Profiling {label}: {cfg['description']} ---")
        info = get_data_info(path)
        data_info[label] = info
        print(f"  Shape     : {info['n_obs']} × {info['n_vars']}")
        print(f"  obs cols  : {info['n_obs_columns']}")
        print(f"  Sparse    : {info.get('sparse', 'N/A')}, dtype={info.get('dtype', 'N/A')}")
        if info.get("sparse"):
            print(f"  Sparsity  : {info.get('sparsity_pct', 'N/A')}%")
        print(f"  Anno cols : {info.get('first_anno_col', 'N/A')}")

    all_results["data_info"] = data_info
    available = [(l, c) for l, c in DATASETS.items() if l in data_info]

    # === 1. scanpy.read_h5ad ===
    print("\n\n" + "=" * 60)
    print("SECTION 1: scanpy.read_h5ad")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]

        def _read(p=path):
            a = sc.read_h5ad(p)
            shape = a.shape
            del a
            gc.collect()
            return shape

        r = run_bench(f"scanpy.read_h5ad [{label}]", _read, iters)
        all_results[f"read_h5ad_{label}"] = r

    # === 2. organize_column_dtypes ===
    print("\n\n" + "=" * 60)
    print("SECTION 2: organize_column_dtypes")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]
        obs = load_obs(path)
        r = run_bench(f"organize_column_dtypes [{label}]", organize_column_dtypes, iters, obs)
        all_results[f"organize_dtypes_{label}"] = r
        del obs
        gc.collect()

    # === 3. get_annotation_columns ===
    print("\n\n" + "=" * 60)
    print("SECTION 3: get_annotation_columns")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]
        obs = load_obs(path)
        r = run_bench(f"get_annotation_columns [{label}]", get_annotation_columns, iters, obs)
        all_results[f"anno_columns_{label}"] = r
        del obs
        gc.collect()

    # === 4. get_pseudotime_columns ===
    print("\n\n" + "=" * 60)
    print("SECTION 4: get_pseudotime_columns")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]
        obs = load_obs(path)
        r = run_bench(f"get_pseudotime_columns [{label}]", get_pseudotime_columns, iters, obs)
        all_results[f"pseudo_columns_{label}"] = r
        del obs
        gc.collect()

    # === 5. Full /columns pipeline ===
    print("\n\n" + "=" * 60)
    print("SECTION 5: /columns Full Pipeline (read + organize + anno + pseudo)")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]
        r = run_bench(f"columns_pipeline [{label}]", columns_pipeline, iters, path)
        all_results[f"columns_pipeline_{label}"] = r

    # === 6. Full /clusters pipeline ===
    print("\n\n" + "=" * 60)
    print("SECTION 6: /clusters Full Pipeline (read + organize + value_counts)")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]
        anno_col = data_info[label].get("first_anno_col")
        if not anno_col:
            print(f"\n  SKIP {label}: no annotation column detected")
            continue
        r = run_bench(f"clusters_pipeline [{label}]", clusters_pipeline, iters, path, anno_col)
        all_results[f"clusters_pipeline_{label}"] = r

    return all_results


# ---------------------------------------------------------------------------
# PART 2: API-level benchmarks
# ---------------------------------------------------------------------------

def run_api_benchmarks():
    all_results = {}

    print("\n\n" + "#" * 60)
    print("# PART 2: API Endpoint Benchmarks (Real scRNA-seq Data)")
    print("#" * 60)

    # Import here to delay app init
    os.environ.pop("TESTING", None)
    from fastapi.testclient import TestClient
    from app.main import app

    print("\n--- Setting up FastAPI TestClient ---")
    client = TestClient(app)

    # Authenticate
    print("--- Authenticating as admin ---")
    response = client.post(
        "/routes/auth/login/access-token",
        data={
            "username": "cellcraft@cellcraft.com",
            "password": "cellcraft2024!",
        },
    )
    if response.status_code != 200:
        print(f"  Login failed ({response.status_code}): {response.text}")
        return {}
    headers = {"Authorization": f"Bearer {response.json()['access_token']}"}
    print("  Auth OK")

    # Only test CRC and GTEx via API (pbmc already measured in previous benchmarks)
    api_datasets = {}
    for label, cfg in DATASETS.items():
        path = cfg["path"]
        if not os.path.isfile(path):
            continue
        # Get annotation column
        adata = sc.read_h5ad(path)
        obs_organized = organize_column_dtypes(adata.obs)
        anno_cols = get_annotation_columns(obs_organized)
        anno_col = anno_cols[0] if anno_cols else None
        del adata
        gc.collect()

        filename = os.path.basename(path)
        api_datasets[label] = {
            "filename": filename,
            "anno_col": anno_col,
            "iterations": cfg["iterations_api"],
        }
        print(f"  {label}: {filename}, anno_col={anno_col}, iters={cfg['iterations_api']}")

    # === /routes/files/columns ===
    print("\n\n" + "=" * 60)
    print("SECTION 1: POST /routes/files/columns")
    print("=" * 60)

    for label, acfg in api_datasets.items():
        iters = acfg["iterations"]
        timer = BenchmarkTimer(f"API /columns [{label}]", iters)

        for i in range(iters):
            start = time.perf_counter()
            resp = client.post(
                "/routes/files/columns",
                json={"file_name": acfg["filename"], "source": "shared"},
                headers=headers,
            )
            elapsed = time.perf_counter() - start
            timer.times.append(elapsed)

            if i == 0:
                if resp.status_code != 200:
                    print(f"\n  ERROR [{label}]: status={resp.status_code}, body={resp.text[:300]}")
                    break
                body = resp.json()
                print(f"\n  [{label}] First response: anno={len(body.get('anno_columns', []))} cols, "
                      f"pseudo={len(body.get('pseudo_columns', []))} cols")

        stats = timer.stats()
        stats["times_ms"] = [round(t * 1000, 2) for t in timer.times]
        if stats.get("iterations", 0) > 1:
            cold = stats["times_ms"][0]
            warm_avg = round(statistics.mean(stats["times_ms"][1:]), 2)
        else:
            cold = stats["times_ms"][0] if stats.get("times_ms") else 0
            warm_avg = cold

        print(f"\n{'='*60}")
        print(f"  {stats['name']}")
        print(f"{'='*60}")
        print(f"  Iterations : {stats['iterations']}")
        print(f"  Mean       : {stats['mean_ms']:.2f} ms")
        print(f"  Median     : {stats['median_ms']:.2f} ms")
        print(f"  Std        : {stats.get('std_s', 0)*1000:.2f} ms")
        print(f"  Min        : {stats['min_s']*1000:.2f} ms")
        print(f"  Max        : {stats['max_s']*1000:.2f} ms")
        print(f"  Cold (1st) : {cold:.2f} ms")
        print(f"  Warm (avg) : {warm_avg:.2f} ms")

        all_results[f"api_columns_{label}"] = stats

    # === /routes/files/clusters ===
    print("\n\n" + "=" * 60)
    print("SECTION 2: POST /routes/files/clusters")
    print("=" * 60)

    for label, acfg in api_datasets.items():
        if not acfg["anno_col"]:
            print(f"\n  SKIP {label}: no annotation column")
            continue

        iters = acfg["iterations"]
        timer = BenchmarkTimer(f"API /clusters [{label}]", iters)

        for i in range(iters):
            start = time.perf_counter()
            resp = client.post(
                "/routes/files/clusters",
                json={
                    "file_name": acfg["filename"],
                    "anno_column": acfg["anno_col"],
                    "source": "shared",
                },
                headers=headers,
            )
            elapsed = time.perf_counter() - start
            timer.times.append(elapsed)

            if i == 0:
                if resp.status_code != 200:
                    print(f"\n  ERROR [{label}]: status={resp.status_code}, body={resp.text[:300]}")
                    break
                body = resp.json()
                print(f"\n  [{label}] First response: {len(body.get('clusters', []))} clusters")

        stats = timer.stats()
        stats["times_ms"] = [round(t * 1000, 2) for t in timer.times]
        if stats.get("iterations", 0) > 1:
            cold = stats["times_ms"][0]
            warm_avg = round(statistics.mean(stats["times_ms"][1:]), 2)
        else:
            cold = stats["times_ms"][0] if stats.get("times_ms") else 0
            warm_avg = cold

        print(f"\n{'='*60}")
        print(f"  {stats['name']}")
        print(f"{'='*60}")
        print(f"  Iterations : {stats['iterations']}")
        print(f"  Mean       : {stats['mean_ms']:.2f} ms")
        print(f"  Median     : {stats['median_ms']:.2f} ms")
        print(f"  Std        : {stats.get('std_s', 0)*1000:.2f} ms")
        print(f"  Min        : {stats['min_s']*1000:.2f} ms")
        print(f"  Max        : {stats['max_s']*1000:.2f} ms")
        print(f"  Cold (1st) : {cold:.2f} ms")
        print(f"  Warm (avg) : {warm_avg:.2f} ms")

        all_results[f"api_clusters_{label}"] = stats

    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "#" * 70)
    print("#  Real scRNA-seq Data Benchmark — Before Redis Caching Baseline")
    print("#" * 70)

    # Part 1
    func_results = run_function_benchmarks()

    # Part 2
    api_results = run_api_benchmarks()

    # Merge
    all_results = {"function_level": func_results, "api_level": api_results}

    # ==================================================================
    # Summary tables
    # ==================================================================
    print("\n\n" + "=" * 70)
    print("FUNCTION-LEVEL SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Benchmark':<50} {'Mean (ms)':>10} {'Median (ms)':>12} {'Peak Mem (MB)':>14}")
    print("-" * 90)
    for key, val in func_results.items():
        if key == "data_info":
            continue
        mem_str = str(val.get("tracemalloc_peak_mb", "N/A"))
        print(f"{val['name']:<50} {val['mean_ms']:>10.2f} {val['median_ms']:>12.2f} {mem_str:>14}")

    print("\n\n" + "=" * 70)
    print("API-LEVEL SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Endpoint':<50} {'Mean (ms)':>10} {'Median (ms)':>12} {'Cold (ms)':>10} {'Warm avg':>10}")
    print("-" * 95)
    for key, val in api_results.items():
        times = val.get("times_ms", [])
        cold = times[0] if times else 0
        warm_avg = round(statistics.mean(times[1:]), 2) if len(times) > 1 else 0
        print(f"{val['name']:<50} {val['mean_ms']:>10.2f} {val['median_ms']:>12.2f} {cold:>10.2f} {warm_avg:>10.2f}")

    # JSON output
    print("\n\n--- RAW JSON ---")
    print(json.dumps(all_results, indent=2))

    return all_results


if __name__ == "__main__":
    main()
