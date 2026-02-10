#!/usr/bin/env python
"""
DataTable benchmarks — Before/After Redis caching.

Measures function-level (load_tab_file) and API-level (POST /datatable/load_data)
performance to establish baselines and validate Redis caching improvements.

Usage (function-level, inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_datatable_benchmarks.py --func

Usage (API-level, inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_datatable_benchmarks.py --api

Usage (both):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_datatable_benchmarks.py
"""
import argparse
import gc
import json
import os
import sys
import time
import statistics
import tracemalloc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATASETS = {
    "pbmc": {
        "path": "/app/tutorials/pbmc_light_1000.h5ad",
        "iterations_func": 5,
        "iterations_api": 10,
        "description": "1K cells, 27 MB",
    },
    "CRC": {
        "path": "/app/tutorials/CRCBraw.h5ad",
        "iterations_func": 3,
        "iterations_api": 5,
        "description": "12.6K cells, 138 MB",
    },
    "GTEx": {
        "path": "/app/tutorials/GTEx_8_tissues.h5ad",
        "iterations_func": 2,
        "iterations_api": 3,
        "description": "209K cells, 1.8 GB",
    },
}

# ---------------------------------------------------------------------------
# Benchmark utilities (reused from run_realdata_benchmarks.py)
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


def print_bench(stats, mem=None):
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


# ---------------------------------------------------------------------------
# PART 1: Function-level benchmarks
# ---------------------------------------------------------------------------

def run_function_benchmarks():
    from app.common.utils.workflow_utils import load_tab_file

    all_results = {}

    print("\n" + "#" * 60)
    print("# PART 1: Function-Level Benchmarks (load_tab_file)")
    print("#" * 60)

    available = []
    for label, cfg in DATASETS.items():
        path = cfg["path"]
        if not os.path.isfile(path):
            print(f"\n  WARNING: {path} not found — skipping {label}")
            continue
        size_mb = round(os.path.getsize(path) / (1024 * 1024), 1)
        print(f"\n--- {label}: {cfg['description']}, file={size_mb} MB ---")
        available.append((label, cfg))

    # === load_tab_file ===
    print("\n\n" + "=" * 60)
    print("SECTION 1: load_tab_file (sc.read_h5ad + obs + UMAP concat)")
    print("=" * 60)

    for label, cfg in available:
        path = cfg["path"]
        iters = cfg["iterations_func"]

        def _load(p=path):
            df = load_tab_file(p)
            shape = df.shape
            del df
            gc.collect()
            return shape

        timer = BenchmarkTimer(f"load_tab_file [{label}]", iters)
        timer.run(_load)
        stats = timer.stats()

        mem_result = measure_memory(_load)
        mem = {k: v for k, v in mem_result.items() if k != "result"}
        df_shape = mem_result["result"]

        print_bench(stats, mem)
        print(f"  DataFrame  : {df_shape[0]} rows × {df_shape[1]} cols")

        all_results[f"load_tab_file_{label}"] = {**stats, **mem, "df_shape": list(df_shape)}

    return all_results


# ---------------------------------------------------------------------------
# PART 2: API-level benchmarks
# ---------------------------------------------------------------------------

def run_api_benchmarks():
    all_results = {}

    print("\n\n" + "#" * 60)
    print("# PART 2: API Endpoint Benchmarks (POST /datatable/load_data)")
    print("#" * 60)

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

    # Filter available datasets
    api_datasets = {}
    for label, cfg in DATASETS.items():
        path = cfg["path"]
        if not os.path.isfile(path):
            continue
        filename = os.path.basename(path)
        api_datasets[label] = {
            "filename": filename,
            "iterations": cfg["iterations_api"],
        }
        print(f"  {label}: {filename}, iters={cfg['iterations_api']}")

    def make_payload(filename, page=1, per_page=20, sort_field=None, sort_type=None, filters=None):
        payload = {
            "file_name": filename,
            "source": "shared",
            "page": page,
            "perPage": per_page,
            "columnFilters": filters or {},
            "sort": [],
        }
        if sort_field and sort_type:
            payload["sort"] = [{"field": sort_field, "type": sort_type}]
        return payload

    # === Scenario 1: Cold / Warm page loads ===
    print("\n\n" + "=" * 60)
    print("SCENARIO 1: Repeated page loads (page 1, same params)")
    print("=" * 60)

    for label, acfg in api_datasets.items():
        iters = acfg["iterations"]
        timer = BenchmarkTimer(f"API load_data repeated [{label}]", iters)
        skip = False

        for i in range(iters):
            payload = make_payload(acfg["filename"])
            start = time.perf_counter()
            try:
                resp = client.post("/routes/datatable/load_data", json=payload, headers=headers)
                elapsed = time.perf_counter() - start
            except Exception as e:
                elapsed = time.perf_counter() - start
                print(f"\n  ERROR [{label}] iter {i}: {type(e).__name__}: {str(e)[:200]}")
                skip = True
                break
            timer.times.append(elapsed)

            if i == 0:
                if resp.status_code != 200:
                    print(f"\n  ERROR [{label}]: status={resp.status_code}, body={resp.text[:300]}")
                    skip = True
                    break
                body = resp.json()
                print(f"\n  [{label}] First response: {len(body.get('rows', []))} rows, "
                      f"{body.get('totalRecords', 0)} total, "
                      f"{len(body.get('columns', []))} columns")

        if skip:
            print(f"  SKIPPING {label} due to error (likely NaN in data)")
            continue

        stats = timer.stats()
        times_ms = [round(t * 1000, 2) for t in timer.times]
        stats["times_ms"] = times_ms
        cold = times_ms[0] if times_ms else 0
        warm_avg = round(statistics.mean(times_ms[1:]), 2) if len(times_ms) > 1 else cold

        print_bench(stats)
        print(f"  Cold (1st) : {cold:.2f} ms")
        print(f"  Warm (avg) : {warm_avg:.2f} ms")

        all_results[f"api_repeated_{label}"] = {**stats, "cold_ms": cold, "warm_avg_ms": warm_avg}

    # === Scenario 2: Page navigation ===
    print("\n\n" + "=" * 60)
    print("SCENARIO 2: Page navigation (page 1 → 2 → 3)")
    print("=" * 60)

    for label, acfg in api_datasets.items():
        times = []
        skip = False
        for page in [1, 2, 3]:
            payload = make_payload(acfg["filename"], page=page)
            start = time.perf_counter()
            try:
                resp = client.post("/routes/datatable/load_data", json=payload, headers=headers)
                elapsed = time.perf_counter() - start
                times.append({"page": page, "ms": round(elapsed * 1000, 2), "status": resp.status_code})
            except Exception as e:
                elapsed = time.perf_counter() - start
                print(f"\n  ERROR [{label}] page {page}: {type(e).__name__}: {str(e)[:200]}")
                skip = True
                break

        if skip:
            print(f"  SKIPPING {label} due to error")
            continue

        print(f"\n  [{label}] Page navigation:")
        for t in times:
            print(f"    Page {t['page']}: {t['ms']:.2f} ms (status={t['status']})")

        all_results[f"api_pages_{label}"] = times

    # === Scenario 3: Sort changes ===
    print("\n\n" + "=" * 60)
    print("SCENARIO 3: Sort changes (same file, different sort)")
    print("=" * 60)

    for label, acfg in api_datasets.items():
        # First get column names from a load
        payload = make_payload(acfg["filename"])
        try:
            resp = client.post("/routes/datatable/load_data", json=payload, headers=headers)
        except Exception as e:
            print(f"\n  SKIPPING {label} sort test: {type(e).__name__}: {str(e)[:200]}")
            continue
        if resp.status_code != 200:
            print(f"\n  SKIPPING {label} sort test: status={resp.status_code}")
            continue
        columns = resp.json().get("columns", [])
        sort_fields = [c["field"] for c in columns[:3]] if columns else []

        times = []
        skip = False
        for field in sort_fields:
            for sort_type in ["asc", "desc"]:
                payload = make_payload(acfg["filename"], sort_field=field, sort_type=sort_type)
                start = time.perf_counter()
                try:
                    resp = client.post("/routes/datatable/load_data", json=payload, headers=headers)
                    elapsed = time.perf_counter() - start
                    times.append({
                        "sort": f"{field}/{sort_type}",
                        "ms": round(elapsed * 1000, 2),
                        "status": resp.status_code,
                    })
                except Exception as e:
                    elapsed = time.perf_counter() - start
                    print(f"\n  ERROR [{label}] sort {field}/{sort_type}: {type(e).__name__}")
                    skip = True
                    break
            if skip:
                break

        if skip:
            print(f"  SKIPPING {label} sort test due to error")
            continue

        print(f"\n  [{label}] Sort changes:")
        for t in times:
            print(f"    {t['sort']}: {t['ms']:.2f} ms (status={t['status']})")

        all_results[f"api_sorts_{label}"] = times

    return all_results


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(func_results, api_results):
    print("\n\n" + "=" * 70)
    print("FUNCTION-LEVEL SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Benchmark':<45} {'Mean (ms)':>10} {'Median (ms)':>12} {'Peak Mem (MB)':>14}")
    print("-" * 85)
    for key, val in func_results.items():
        mem_str = str(val.get("tracemalloc_peak_mb", "N/A"))
        print(f"{val['name']:<45} {val['mean_ms']:>10.2f} {val['median_ms']:>12.2f} {mem_str:>14}")

    print("\n\n" + "=" * 70)
    print("API-LEVEL SUMMARY TABLE")
    print("=" * 70)

    # Repeated loads summary
    print(f"\n{'Endpoint':<45} {'Mean (ms)':>10} {'Cold (ms)':>10} {'Warm avg':>10}")
    print("-" * 80)
    for key, val in api_results.items():
        if key.startswith("api_repeated_"):
            print(f"{val['name']:<45} {val['mean_ms']:>10.2f} {val.get('cold_ms', 0):>10.2f} {val.get('warm_avg_ms', 0):>10.2f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="DataTable benchmark suite")
    parser.add_argument("--func", action="store_true", help="Run function-level benchmarks only")
    parser.add_argument("--api", action="store_true", help="Run API-level benchmarks only")
    args = parser.parse_args()

    run_func = args.func or (not args.func and not args.api)
    run_api = args.api or (not args.func and not args.api)

    print("\n" + "#" * 70)
    print("#  DataTable Benchmark — P2 Redis Caching")
    print("#" * 70)

    func_results = {}
    api_results = {}

    if run_func:
        func_results = run_function_benchmarks()

    if run_api:
        api_results = run_api_benchmarks()

    # Summary
    print_summary(func_results, api_results)

    # Combined JSON output
    all_results = {"function_level": func_results, "api_level": api_results}
    print("\n\n--- RAW JSON ---")
    print(json.dumps(all_results, indent=2, default=str))

    return all_results


if __name__ == "__main__":
    main()
