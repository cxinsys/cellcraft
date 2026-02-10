"""API-level test for h5ad metadata caching.

Tests cache by making repeated API calls through the Uvicorn server.
The in-process cache lives in the Uvicorn worker, so we can't clear it
from this script. Instead we rely on:
  - 1st call: cache miss (server just started or file not yet cached)
  - 2nd+ calls: cache hit (in-process cache populated by 1st call)

To test cold performance, restart the backend container first:
  docker restart cellcraft-backend-1 && sleep 3
"""
import httpx
import time
import os
import sys

BASE = "http://localhost:8000/routes"


def login():
    resp = httpx.post(
        f"{BASE}/auth/login/access-token",
        data={"username": "cellcraft@cellcraft.com", "password": "cellcraft2024!"},
    )
    assert resp.status_code == 200, f"Login failed: {resp.status_code}"
    return resp.json()["access_token"]


def timed_post(url, json, headers, timeout=30):
    t0 = time.perf_counter()
    r = httpx.post(url, json=json, headers=headers, timeout=timeout)
    elapsed = (time.perf_counter() - t0) * 1000
    return r, elapsed


def main():
    token = login()
    headers = {"Authorization": f"Bearer {token}"}
    print("Login OK\n")

    datasets = [
        ("pbmc_light_1000.h5ad", "./tutorials/pbmc_light_1000.h5ad"),
        ("CRCBraw.h5ad", "./tutorials/CRCBraw.h5ad"),
    ]

    for fname, fpath in datasets:
        if not os.path.isfile(fpath):
            print(f"[SKIP] {fname} not found at {fpath}")
            continue

        print(f"=== {fname} ===")
        payload = {"file_name": fname, "source": "shared"}

        # 1st call (may be cold or warm depending on server state)
        r1, t1 = timed_post(f"{BASE}/files/columns", payload, headers, timeout=60)
        assert r1.status_code == 200, f"/columns failed: {r1.status_code} {r1.text}"

        # 2nd call (should be warm/cached)
        r2, t2 = timed_post(f"{BASE}/files/columns", payload, headers, timeout=60)
        assert r2.status_code == 200
        assert r1.json() == r2.json(), "Results differ between calls!"

        # 3rd call to confirm consistency
        r3, t3 = timed_post(f"{BASE}/files/columns", payload, headers, timeout=60)
        assert r3.status_code == 200
        assert r1.json() == r3.json()

        print(f"  /columns: call1={t1:.0f}ms  call2={t2:.0f}ms  call3={t3:.0f}ms")

        # Clusters
        anno_col = r1.json()["anno_columns"][0]
        payload_cl = {"file_name": fname, "anno_column": anno_col, "source": "shared"}

        rc1, tc1 = timed_post(f"{BASE}/files/clusters", payload_cl, headers, timeout=60)
        rc2, tc2 = timed_post(f"{BASE}/files/clusters", payload_cl, headers, timeout=60)
        assert rc1.status_code == 200 and rc2.status_code == 200
        assert rc1.json() == rc2.json()

        print(f"  /clusters ({anno_col}): call1={tc1:.0f}ms  call2={tc2:.0f}ms")
        print()

    print("All API cache tests PASSED!")


if __name__ == "__main__":
    main()
