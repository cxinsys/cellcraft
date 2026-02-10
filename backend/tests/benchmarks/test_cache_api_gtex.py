"""API-level cache benchmark for all 3 datasets including GTEx.

Run after: docker restart cellcraft-backend-1 && sleep 5
"""
import httpx
import time
import os

BASE = "http://localhost:8000/routes"


def login():
    resp = httpx.post(
        f"{BASE}/auth/login/access-token",
        data={"username": "cellcraft@cellcraft.com", "password": "cellcraft2024!"},
    )
    assert resp.status_code == 200, f"Login failed: {resp.status_code}"
    return resp.json()["access_token"]


def timed_post(url, json, headers, timeout=120):
    t0 = time.perf_counter()
    r = httpx.post(url, json=json, headers=headers, timeout=timeout)
    elapsed = (time.perf_counter() - t0) * 1000
    return r, elapsed


def main():
    token = login()
    headers = {"Authorization": f"Bearer {token}"}
    print("Login OK\n")

    datasets = [
        ("pbmc_light_1000.h5ad", "./tutorials/pbmc_light_1000.h5ad", "pbmc (1K cells, 27MB)"),
        ("CRCBraw.h5ad", "./tutorials/CRCBraw.h5ad", "CRC (12.6K cells, 138MB)"),
        ("GTEx_8_tissues.h5ad", "./tutorials/GTEx_8_tissues.h5ad", "GTEx (209K cells, 1.8GB)"),
    ]

    print(f"{'Dataset':<30} {'Endpoint':<12} {'Cold (ms)':>10} {'Warm1 (ms)':>11} {'Warm2 (ms)':>11} {'Speedup':>8}")
    print("-" * 90)

    for fname, fpath, label in datasets:
        if not os.path.isfile(fpath):
            print(f"[SKIP] {label} — file not found")
            continue

        payload = {"file_name": fname, "source": "shared"}

        # /columns: cold
        r1, t_cold = timed_post(f"{BASE}/files/columns", payload, headers)
        assert r1.status_code == 200, f"/columns failed: {r1.status_code} {r1.text}"

        # /columns: warm x2
        r2, t_warm1 = timed_post(f"{BASE}/files/columns", payload, headers)
        r3, t_warm2 = timed_post(f"{BASE}/files/columns", payload, headers)
        assert r2.status_code == 200 and r3.status_code == 200
        assert r1.json() == r2.json() == r3.json()

        speedup = t_cold / t_warm1
        print(f"{label:<30} {'/ columns':<12} {t_cold:>10.0f} {t_warm1:>11.0f} {t_warm2:>11.0f} {speedup:>7.1f}x")

        # /clusters
        anno_col = r1.json()["anno_columns"][0]
        payload_cl = {"file_name": fname, "anno_column": anno_col, "source": "shared"}

        rc1, tc_cold = timed_post(f"{BASE}/files/clusters", payload_cl, headers)
        rc2, tc_warm1 = timed_post(f"{BASE}/files/clusters", payload_cl, headers)
        rc3, tc_warm2 = timed_post(f"{BASE}/files/clusters", payload_cl, headers)
        assert rc1.status_code == 200 and rc2.status_code == 200 and rc3.status_code == 200
        assert rc1.json() == rc2.json() == rc3.json()

        speedup_cl = tc_cold / tc_warm1
        print(f"{'':30} {'/ clusters':<12} {tc_cold:>10.0f} {tc_warm1:>11.0f} {tc_warm2:>11.0f} {speedup_cl:>7.1f}x")
        print()

    print("All tests PASSED!")


if __name__ == "__main__":
    main()
