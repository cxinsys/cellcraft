#!/usr/bin/env python
"""
Plugin List benchmarks — Before/After Redis caching + code optimizations.

Measures GET /plugin/list performance across three scenarios:
  1. Before (original code path simulation): cache disabled, old loop pattern
  2. Cold  (optimized, cache miss): cache cleared before each call
  3. Warm  (cache hit): cached response served from Redis

Usage (inside backend container):
    docker exec cellcraft-backend-1 python tests/benchmarks/run_plugin_list_benchmarks.py
"""
import json
import os
import sys
import time
import statistics

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


# ---------------------------------------------------------------------------
# Benchmark utilities
# ---------------------------------------------------------------------------

class BenchmarkTimer:
    def __init__(self, name: str):
        self.name = name
        self.times = []

    def record(self, elapsed: float):
        self.times.append(elapsed)

    def stats(self):
        if not self.times:
            return {"name": self.name, "iterations": 0}
        s = sorted(self.times)
        n = len(s)
        return {
            "name": self.name,
            "iterations": n,
            "mean_ms": round(statistics.mean(s) * 1000, 2),
            "median_ms": round(statistics.median(s) * 1000, 2),
            "std_ms": round(statistics.stdev(s) * 1000, 2) if n > 1 else 0.0,
            "min_ms": round(s[0] * 1000, 2),
            "max_ms": round(s[-1] * 1000, 2),
            "times_ms": [round(t * 1000, 2) for t in self.times],
        }


def print_bench(stats):
    print(f"\n{'='*60}")
    print(f"  {stats['name']}")
    print(f"{'='*60}")
    print(f"  Iterations : {stats['iterations']}")
    print(f"  Mean       : {stats['mean_ms']:.2f} ms")
    print(f"  Median     : {stats['median_ms']:.2f} ms")
    print(f"  Std        : {stats.get('std_ms', 0):.2f} ms")
    print(f"  Min        : {stats['min_ms']:.2f} ms")
    print(f"  Max        : {stats['max_ms']:.2f} ms")
    if stats.get("times_ms"):
        print(f"  All times  : {stats['times_ms']}")


# ---------------------------------------------------------------------------
# Setup: FastAPI TestClient + Auth
# ---------------------------------------------------------------------------

def setup_client():
    os.environ.pop("TESTING", None)
    from fastapi.testclient import TestClient
    from app.main import app

    print("\n--- Setting up FastAPI TestClient ---")
    client = TestClient(app)

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
        sys.exit(1)

    token = response.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}
    print("  Auth OK")

    # Extract user_id from token for cache operations
    from app.auth.deps import get_current_active_user
    import jwt
    from app.core.config import settings
    payload = jwt.decode(token, settings.SECRET_KEY, algorithms=["HS256"])
    user_id = payload.get("sub")
    print(f"  User ID: {user_id}")

    return client, headers, user_id


# ---------------------------------------------------------------------------
# Redis cache helpers
# ---------------------------------------------------------------------------

def clear_plugin_cache():
    """Clear all plugin list cache entries from Redis."""
    try:
        from app.shared.redis import cache_delete_pattern
        count = cache_delete_pattern("plugin:list:*")
        return count
    except Exception as e:
        print(f"  WARNING: Failed to clear plugin cache: {e}")
        return 0


def check_cache_exists(user_id):
    """Check if cache exists for user."""
    try:
        from app.shared.redis import cache_get_bytes
        data = cache_get_bytes(f"plugin:list:{user_id}")
        return data is not None
    except Exception:
        return False


def get_cache_size(user_id):
    """Get size of cached data in bytes."""
    try:
        from app.shared.redis import cache_get_bytes
        data = cache_get_bytes(f"plugin:list:{user_id}")
        return len(data) if data else 0
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# Scenario 1: "Before" — Original code path simulation
# ---------------------------------------------------------------------------

def run_before_benchmark(client, headers, iterations=5):
    """
    Simulate original code path by clearing cache before every call.

    Note: The code optimizations (get_user_task outside loop, Docker client reuse)
    are already in place. This measures the "cold" path with optimizations applied.
    For true "before" comparison, we also run the old-style function directly.
    """
    print("\n\n" + "#" * 60)
    print("# SCENARIO 1: Before (cache disabled — clear before each call)")
    print("# This simulates pre-cache behavior. Code optimizations are")
    print("# already applied, so this is optimized-cold, not true-original.")
    print("#" * 60)

    timer = BenchmarkTimer("GET /plugin/list [cache disabled]")

    for i in range(iterations):
        clear_plugin_cache()

        start = time.perf_counter()
        resp = client.get("/routes/plugin/list", headers=headers)
        elapsed = time.perf_counter() - start
        timer.record(elapsed)

        if i == 0:
            if resp.status_code != 200:
                print(f"  ERROR: status={resp.status_code}, body={resp.text[:300]}")
                return {}
            body = resp.json()
            plugin_count = len(body.get("plugins", []))
            print(f"\n  Response: {plugin_count} plugins, status={resp.status_code}")

        # Clear cache after each call to force next call to be cold
        clear_plugin_cache()

    stats = timer.stats()
    print_bench(stats)
    return stats


# ---------------------------------------------------------------------------
# Scenario 1b: True "Before" — Old code path with N queries + N docker calls
# ---------------------------------------------------------------------------

def run_true_before_benchmark(db_session, user_id, iterations=5):
    """
    Execute the original list_plugins logic (before code optimizations):
    - get_user_task() called N times inside loop
    - docker.from_env() called N times inside loop

    This gives us the true "before" baseline.
    """
    print("\n\n" + "#" * 60)
    print("# SCENARIO 1b: True Before (original code path, no optimizations)")
    print("# Reproduces N×get_user_task + N×docker.from_env loop pattern")
    print("#" * 60)

    from app.plugin import crud as crud_plugin
    from app.task.crud import get_user_task
    from celery.result import AsyncResult
    from datetime import datetime
    import docker

    timer = BenchmarkTimer("list_plugins [original code path]")

    for iteration in range(iterations):
        start = time.perf_counter()

        # === Original code path (pre-optimization) ===
        plugins = crud_plugin.get_plugins(db_session)
        plugin_list = []

        for plugin in plugins:
            plugin_dict = plugin.__dict__.copy()
            plugin_dict['users'] = [user.__dict__ for user in plugin.users]

            if 'rules' in plugin_dict and plugin_dict['rules']:
                rules_dict = plugin_dict['rules']
                rules_array = [rules_dict[str(i)] for i in range(len(rules_dict))]
                plugin_dict['rules'] = rules_array

            # Original: check_plugin_docker_image per plugin (creates client each time)
            try:
                if plugin.source == "local":
                    docker_client = None
                    try:
                        docker_client = docker.from_env()
                        image_tag = f"plugin-{plugin.name.lower()}"
                        docker_client.images.get(image_tag)
                        plugin_dict['docker_image_exists'] = True
                    except docker.errors.ImageNotFound:
                        plugin_dict['docker_image_exists'] = False
                    except Exception:
                        plugin_dict['docker_image_exists'] = False
                    finally:
                        if docker_client:
                            docker_client.close()
                else:
                    plugin_dict['docker_image_exists'] = False
            except Exception:
                plugin_dict['docker_image_exists'] = False

            # Original: get_user_task called EVERY iteration (N times)
            try:
                user_tasks = get_user_task(db_session, int(user_id))
                plugin_build_tasks = [
                    task for task in user_tasks
                    if task.task_type == "plugin_build" and task.plugin_name == plugin.name
                ]
                if plugin_build_tasks:
                    latest_task = max(plugin_build_tasks,
                                     key=lambda x: x.start_time or datetime.min.replace(tzinfo=None))
                    task_result = AsyncResult(latest_task.task_id)
                    plugin_dict['latest_build'] = {
                        "task_id": latest_task.task_id,
                        "status": task_result.state,
                        "start_time": latest_task.start_time.isoformat() if latest_task.start_time else None,
                        "end_time": latest_task.end_time.isoformat() if latest_task.end_time else None,
                    }
                else:
                    plugin_dict['latest_build'] = None
            except Exception:
                plugin_dict['latest_build'] = None

            if plugin.source == "official":
                plugin_dict['current_version'] = plugin.version or "1.0"
                plugin_dict['available_versions'] = []
            else:
                plugin_dict['current_version'] = plugin.version or "local"
                plugin_dict['available_versions'] = []

            plugin_list.append(plugin_dict)

        response = {"plugins": plugin_list}
        elapsed = time.perf_counter() - start
        timer.record(elapsed)

        if iteration == 0:
            print(f"\n  Response: {len(plugin_list)} plugins")

    stats = timer.stats()
    print_bench(stats)
    return stats


# ---------------------------------------------------------------------------
# Scenario 2: Cold (optimized, cache miss)
# ---------------------------------------------------------------------------

def run_cold_benchmark(client, headers, user_id, iterations=5):
    """Measure optimized code with cache miss (clear before each call)."""
    print("\n\n" + "#" * 60)
    print("# SCENARIO 2: Cold (optimized code, cache miss)")
    print("# Cache is cleared before each call")
    print("#" * 60)

    timer = BenchmarkTimer("GET /plugin/list [cold — cache miss]")

    for i in range(iterations):
        clear_plugin_cache()

        start = time.perf_counter()
        resp = client.get("/routes/plugin/list", headers=headers)
        elapsed = time.perf_counter() - start
        timer.record(elapsed)

        if i == 0 and resp.status_code != 200:
            print(f"  ERROR: status={resp.status_code}")
            return {}

        # Verify cache was created
        if i == 0:
            exists = check_cache_exists(user_id)
            size = get_cache_size(user_id)
            print(f"\n  Cache created: {exists}, size: {size} bytes ({size/1024:.1f} KB)")

        clear_plugin_cache()

    stats = timer.stats()
    print_bench(stats)
    return stats


# ---------------------------------------------------------------------------
# Scenario 3: Warm (cache hit)
# ---------------------------------------------------------------------------

def run_warm_benchmark(client, headers, user_id, iterations=10):
    """Measure cache hit performance."""
    print("\n\n" + "#" * 60)
    print("# SCENARIO 3: Warm (cache hit)")
    print("# First call populates cache, subsequent calls are cache hits")
    print("#" * 60)

    # Ensure cache is populated
    clear_plugin_cache()
    resp = client.get("/routes/plugin/list", headers=headers)
    if resp.status_code != 200:
        print(f"  ERROR populating cache: status={resp.status_code}")
        return {}

    exists = check_cache_exists(user_id)
    size = get_cache_size(user_id)
    print(f"\n  Cache populated: {exists}, size: {size} bytes ({size/1024:.1f} KB)")

    timer = BenchmarkTimer("GET /plugin/list [warm — cache hit]")

    for i in range(iterations):
        start = time.perf_counter()
        resp = client.get("/routes/plugin/list", headers=headers)
        elapsed = time.perf_counter() - start
        timer.record(elapsed)

    stats = timer.stats()
    print_bench(stats)
    return stats


# ---------------------------------------------------------------------------
# Scenario 4: Cache invalidation + re-population
# ---------------------------------------------------------------------------

def run_invalidation_benchmark(client, headers, user_id, iterations=5):
    """Measure invalidation + re-population cycle."""
    print("\n\n" + "#" * 60)
    print("# SCENARIO 4: Invalidation + Re-population cycle")
    print("# Simulates associate/dissociate → next list call")
    print("#" * 60)

    timer_invalidate = BenchmarkTimer("invalidate_all_plugin_cache()")
    timer_repopulate = BenchmarkTimer("GET /plugin/list [after invalidation]")

    # Ensure cache exists
    clear_plugin_cache()
    client.get("/routes/plugin/list", headers=headers)

    for i in range(iterations):
        # Measure invalidation
        start = time.perf_counter()
        from app.plugin.cache import invalidate_all_plugin_cache
        invalidate_all_plugin_cache()
        elapsed = time.perf_counter() - start
        timer_invalidate.record(elapsed)

        # Measure re-population
        start = time.perf_counter()
        resp = client.get("/routes/plugin/list", headers=headers)
        elapsed = time.perf_counter() - start
        timer_repopulate.record(elapsed)

    stats_inv = timer_invalidate.stats()
    stats_rep = timer_repopulate.stats()
    print_bench(stats_inv)
    print_bench(stats_rep)
    return {"invalidation": stats_inv, "repopulation": stats_rep}


# ---------------------------------------------------------------------------
# Scenario 5: Redis fallback (graceful degradation)
# ---------------------------------------------------------------------------

def run_fallback_check(client, headers):
    """Verify endpoint works when Redis is unreachable (manual check instruction)."""
    print("\n\n" + "#" * 60)
    print("# SCENARIO 5: Redis fallback verification")
    print("#" * 60)
    print("\n  To verify Redis fallback, run manually:")
    print("    1. docker stop cellcraft-redis-1")
    print("    2. curl -H 'Authorization: Bearer <token>' http://localhost:8000/routes/plugin/list")
    print("    3. Verify response is 200 with plugins data")
    print("    4. docker start cellcraft-redis-1")
    print("\n  (Automated test skipped — stopping Redis affects other services)")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(results):
    print("\n\n" + "=" * 70)
    print("SUMMARY — P3: Plugin List Cache Performance")
    print("=" * 70)

    print(f"\n{'Scenario':<50} {'Mean (ms)':>10} {'Median (ms)':>12}")
    print("-" * 75)

    for key in ["true_before", "cold", "warm"]:
        if key in results and results[key].get("iterations", 0) > 0:
            r = results[key]
            print(f"{r['name']:<50} {r['mean_ms']:>10.2f} {r['median_ms']:>12.2f}")

    # Speedup calculations
    print("\n" + "-" * 75)
    before = results.get("true_before", {})
    cold = results.get("cold", {})
    warm = results.get("warm", {})

    if before.get("mean_ms") and cold.get("mean_ms"):
        speedup = before["mean_ms"] / cold["mean_ms"]
        print(f"  Code optimization speedup (before→cold): {speedup:.1f}x")
        print(f"    ({before['mean_ms']:.0f}ms → {cold['mean_ms']:.0f}ms)")

    if before.get("mean_ms") and warm.get("mean_ms"):
        speedup = before["mean_ms"] / warm["mean_ms"]
        print(f"  Cache hit speedup (before→warm):          {speedup:.0f}x")
        print(f"    ({before['mean_ms']:.0f}ms → {warm['mean_ms']:.0f}ms)")

    if cold.get("mean_ms") and warm.get("mean_ms"):
        speedup = cold["mean_ms"] / warm["mean_ms"]
        print(f"  Cache hit speedup (cold→warm):            {speedup:.0f}x")
        print(f"    ({cold['mean_ms']:.0f}ms → {warm['mean_ms']:.0f}ms)")

    # Invalidation info
    inv = results.get("invalidation", {})
    if isinstance(inv, dict) and "invalidation" in inv:
        inv_stats = inv["invalidation"]
        rep_stats = inv["repopulation"]
        print(f"\n  Invalidation time: {inv_stats.get('mean_ms', 0):.2f}ms avg")
        print(f"  Re-population time: {rep_stats.get('mean_ms', 0):.2f}ms avg")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "#" * 70)
    print("#  P3: Plugin List Cache Benchmark")
    print("#  Before/After Redis Caching + Code Optimizations")
    print("#" * 70)

    client, headers, user_id = setup_client()

    results = {}

    # Scenario 1b: True before (original code path)
    try:
        from app.db.session import get_db_session
        with get_db_session() as db:
            results["true_before"] = run_true_before_benchmark(db, user_id, iterations=5)
    except Exception as e:
        print(f"\n  WARNING: True before benchmark failed: {e}")
        import traceback
        traceback.print_exc()

    # Scenario 2: Cold (optimized, cache miss)
    results["cold"] = run_cold_benchmark(client, headers, user_id, iterations=5)

    # Scenario 3: Warm (cache hit)
    results["warm"] = run_warm_benchmark(client, headers, user_id, iterations=10)

    # Scenario 4: Invalidation cycle
    results["invalidation"] = run_invalidation_benchmark(client, headers, user_id, iterations=5)

    # Scenario 5: Fallback info
    run_fallback_check(client, headers)

    # Summary
    print_summary(results)

    # Clean up
    clear_plugin_cache()

    # Raw JSON
    print("\n\n--- RAW JSON ---")
    print(json.dumps(results, indent=2, default=str))

    return results


if __name__ == "__main__":
    main()
