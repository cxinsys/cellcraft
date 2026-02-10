#!/usr/bin/env python
"""
H5AD API endpoint benchmarks — Before Redis caching baseline.

Standalone script (no pytest) to measure actual HTTP response times
through FastAPI TestClient for /routes/files/columns and /routes/files/clusters.

Usage (inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_api_benchmarks.py
"""
import gc
import json
import os
import sys
import time
import statistics
import tempfile

import numpy as np
import pandas as pd
import anndata as ad

# Ensure TESTING is NOT set so the app connects to the real DB (db:5432)
os.environ.pop("TESTING", None)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from fastapi.testclient import TestClient
from app.main import app

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ITERATIONS = 10
REAL_H5AD_NAME = "pbmc_light_1000.h5ad"
REAL_H5AD_PATH = "/app/tutorials/pbmc_light_1000.h5ad"

# ---------------------------------------------------------------------------
# Benchmark timer (inline to avoid import issues)
# ---------------------------------------------------------------------------

class BenchmarkTimer:
    def __init__(self, name, iterations=ITERATIONS):
        self.name = name
        self.iterations = iterations
        self.times = []

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
            "mean_ms": round(statistics.mean(s) * 1000, 2),
            "median_ms": round(statistics.median(s) * 1000, 2),
            "cold_ms": round(s[0] * 1000, 2) if n > 0 else 0,
        }


def print_stats(stats):
    print(f"\n{'='*60}")
    print(f"  {stats['name']}")
    print(f"{'='*60}")
    print(f"  Iterations : {stats['iterations']}")
    print(f"  Mean       : {stats['mean_ms']:.2f} ms")
    print(f"  Median     : {stats['median_ms']:.2f} ms")
    print(f"  Std        : {stats.get('std_s', 0)*1000:.2f} ms")
    print(f"  Min        : {stats['min_s']*1000:.2f} ms")
    print(f"  Max        : {stats['max_s']*1000:.2f} ms")
    if stats['iterations'] > 1:
        cold = stats['times_ms'][0]
        warm_avg = statistics.mean(stats['times_ms'][1:])
        print(f"  Cold (1st) : {cold:.2f} ms")
        print(f"  Warm (avg) : {warm_avg:.2f} ms")


# ---------------------------------------------------------------------------
# Prepare test data (small & large mock h5ad as user-owned files)
# ---------------------------------------------------------------------------

def create_test_h5ad(n_obs, n_vars, path):
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


# ---------------------------------------------------------------------------
# Authentication helper
# ---------------------------------------------------------------------------

def get_auth_token(client):
    """Login as admin user and return auth headers."""
    response = client.post(
        "/routes/auth/login/access-token",
        data={
            "username": "cellcraft@cellcraft.com",
            "password": "cellcraft2024!",
        },
    )
    if response.status_code != 200:
        print(f"  Login failed ({response.status_code}): {response.text}")
        print("  Trying to create admin user or use existing credentials...")
        return None
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------

def main():
    all_results = {}

    print("\n" + "#" * 60)
    print("# H5AD API Endpoint Benchmark — Before Redis Caching Baseline")
    print("#" * 60)

    # --- Setup TestClient ---
    print("\n--- Setting up FastAPI TestClient ---")
    client = TestClient(app)

    # --- Authenticate ---
    print("--- Authenticating as admin ---")
    headers = get_auth_token(client)
    if not headers:
        print("ERROR: Cannot authenticate. Aborting.")
        return

    print("  Auth OK")

    # --- Prepare test data ---
    # For small/large: create h5ad in tutorials/ dir and use source=shared
    # This avoids needing file DB records
    tmpdir = tempfile.mkdtemp(prefix="bench_api_", dir="/app/tutorials")
    small_name = "bench_small.h5ad"
    large_name = "bench_large.h5ad"
    small_path = os.path.join("/app/tutorials", small_name)
    large_path = os.path.join("/app/tutorials", large_name)

    print("\n--- Preparing test data ---")
    small_size = create_test_h5ad(100, 50, small_path)
    print(f"  Small (100×50)   : {small_size:.2f} MB")
    large_size = create_test_h5ad(10_000, 2_000, large_path)
    print(f"  Large (10k×2k)   : {large_size:.2f} MB")
    real_exists = os.path.isfile(REAL_H5AD_PATH)
    if real_exists:
        real_size = os.path.getsize(REAL_H5AD_PATH) / (1024 * 1024)
        print(f"  Real (pbmc 1k)   : {real_size:.2f} MB")
    else:
        print("  Real (pbmc 1k)   : NOT FOUND — skipping")

    # Define test datasets: (label, filename, anno_column_for_clusters)
    datasets = [
        ("small", small_name, "cluster"),
        ("large", large_name, "cluster"),
    ]
    if real_exists:
        # Discover annotation column from real data first
        import scanpy as sc
        from app.common.utils.h5ad_utils import organize_column_dtypes, get_annotation_columns
        adata = sc.read_h5ad(REAL_H5AD_PATH)
        obs = organize_column_dtypes(adata.obs)
        anno_cols = get_annotation_columns(obs)
        real_anno_col = anno_cols[0] if anno_cols else "cluster"
        del adata
        gc.collect()
        print(f"  Real anno column : {real_anno_col}")
        datasets.append(("real", REAL_H5AD_NAME, real_anno_col))

    # ==================================================================
    # /routes/files/columns benchmarks
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 1: POST /routes/files/columns")
    print("=" * 60)

    for label, filename, _ in datasets:
        timer = BenchmarkTimer(f"API /columns [{label}]", ITERATIONS)

        for i in range(ITERATIONS):
            start = time.perf_counter()
            resp = client.post(
                "/routes/files/columns",
                json={"file_name": filename, "source": "shared"},
                headers=headers,
            )
            elapsed = time.perf_counter() - start
            timer.times.append(elapsed)

            if i == 0:
                # Verify first response
                if resp.status_code != 200:
                    print(f"\n  ERROR [{label}]: status={resp.status_code}, body={resp.text[:200]}")
                    break
                body = resp.json()
                print(f"\n  [{label}] First response: anno={len(body.get('anno_columns', []))} cols, "
                      f"pseudo={len(body.get('pseudo_columns', []))} cols")

        stats = timer.stats()
        stats["times_ms"] = [round(t * 1000, 2) for t in timer.times]
        print_stats(stats)
        all_results[f"api_columns_{label}"] = stats

    # ==================================================================
    # /routes/files/clusters benchmarks
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SECTION 2: POST /routes/files/clusters")
    print("=" * 60)

    for label, filename, anno_col in datasets:
        timer = BenchmarkTimer(f"API /clusters [{label}]", ITERATIONS)

        for i in range(ITERATIONS):
            start = time.perf_counter()
            resp = client.post(
                "/routes/files/clusters",
                json={"file_name": filename, "anno_column": anno_col, "source": "shared"},
                headers=headers,
            )
            elapsed = time.perf_counter() - start
            timer.times.append(elapsed)

            if i == 0:
                if resp.status_code != 200:
                    print(f"\n  ERROR [{label}]: status={resp.status_code}, body={resp.text[:200]}")
                    break
                body = resp.json()
                print(f"\n  [{label}] First response: {len(body.get('clusters', []))} clusters")

        stats = timer.stats()
        stats["times_ms"] = [round(t * 1000, 2) for t in timer.times]
        print_stats(stats)
        all_results[f"api_clusters_{label}"] = stats

    # ==================================================================
    # Summary
    # ==================================================================
    print("\n\n" + "=" * 60)
    print("SUMMARY TABLE")
    print("=" * 60)
    print(f"{'Endpoint':<40} {'Mean (ms)':>10} {'Median (ms)':>12} {'Cold (ms)':>10} {'Warm avg (ms)':>14}")
    print("-" * 90)
    for key, val in all_results.items():
        times = val.get("times_ms", [])
        cold = times[0] if times else 0
        warm_avg = round(statistics.mean(times[1:]), 2) if len(times) > 1 else 0
        print(f"{val['name']:<40} {val['mean_ms']:>10.2f} {val['median_ms']:>12.2f} {cold:>10.2f} {warm_avg:>14.2f}")

    # JSON output
    print("\n--- RAW JSON ---")
    print(json.dumps(all_results, indent=2))

    # Cleanup temp files
    try:
        os.unlink(small_path)
        os.unlink(large_path)
        os.rmdir(tmpdir)
    except Exception:
        pass

    return all_results


if __name__ == "__main__":
    main()
