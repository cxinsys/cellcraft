"""
H5AD function-level benchmarks — Before Redis caching baseline.

Measures scanpy.read_h5ad, organize_column_dtypes, get_annotation_columns,
get_pseudotime_columns, and full /columns + /clusters pipeline simulation.

Run inside backend container:
    docker exec cellcraft-backend-1 \
        python -m pytest tests/benchmarks/test_h5ad_function_benchmarks.py \
        -v -s -m benchmark --no-header --override-ini="addopts=" 2>&1
"""
import gc
import json
import os
import tempfile
import time

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from app.datatable.h5ad import (
    get_annotation_columns,
    get_pseudotime_columns,
    organize_column_dtypes,
)

from .utils import BenchmarkTimer, MemoryProfiler, format_benchmark_result, get_real_h5ad_path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ITERATIONS = 10
REAL_H5AD = get_real_h5ad_path()

# ---------------------------------------------------------------------------
# Fixtures — create h5ad files once per module for speed
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def small_h5ad(tmp_path_factory):
    """100 cells × 50 genes"""
    path = str(tmp_path_factory.mktemp("small") / "small.h5ad")
    n_obs, n_vars = 100, 50
    adata = ad.AnnData(
        X=np.random.rand(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": np.random.choice(["TypeA", "TypeB", "TypeC"], n_obs),
            "cluster": np.random.choice(["C1", "C2"], n_obs),
            "pseudotime": np.random.rand(n_obs),
        }),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )
    adata.obsm["X_umap"] = np.random.rand(n_obs, 2)
    adata.write_h5ad(path)
    del adata
    gc.collect()
    return path


@pytest.fixture(scope="module")
def large_h5ad(tmp_path_factory):
    """10 000 cells × 2 000 genes"""
    path = str(tmp_path_factory.mktemp("large") / "large.h5ad")
    n_obs, n_vars = 10_000, 2_000
    adata = ad.AnnData(
        X=np.random.rand(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": np.random.choice(["TypeA", "TypeB", "TypeC", "TypeD"], n_obs),
            "cluster": np.random.choice([f"C{i}" for i in range(10)], n_obs),
            "pseudotime": np.random.rand(n_obs),
        }),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )
    adata.obsm["X_umap"] = np.random.rand(n_obs, 2)
    adata.write_h5ad(path)
    del adata
    gc.collect()
    return path


# ===================================================================
# scanpy.read_h5ad benchmarks
# ===================================================================

@pytest.mark.benchmark
class TestBenchScanpyReadH5ad:

    def test_bench_scanpy_read_h5ad_small(self, small_h5ad):
        """scanpy.read_h5ad — 100×50"""
        timer = BenchmarkTimer("scanpy.read_h5ad [small 100×50]", ITERATIONS)

        def _read():
            adata = sc.read_h5ad(small_h5ad)
            result = adata.shape
            del adata
            gc.collect()
            return result

        timer.run(_read)
        mem = MemoryProfiler.measure(lambda: sc.read_h5ad(small_h5ad))
        stats = timer.stats()

        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

        assert stats["mean_s"] < 30, "read_h5ad small should complete within 30s"

    def test_bench_scanpy_read_h5ad_large(self, large_h5ad):
        """scanpy.read_h5ad — 10000×2000"""
        timer = BenchmarkTimer("scanpy.read_h5ad [large 10k×2k]", ITERATIONS)

        def _read():
            adata = sc.read_h5ad(large_h5ad)
            result = adata.shape
            del adata
            gc.collect()
            return result

        timer.run(_read)
        mem = MemoryProfiler.measure(lambda: sc.read_h5ad(large_h5ad))
        stats = timer.stats()

        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

        assert stats["mean_s"] < 60, "read_h5ad large should complete within 60s"

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_scanpy_read_h5ad_real(self):
        """scanpy.read_h5ad — pbmc_light_1000.h5ad (~28.7 MB)"""
        timer = BenchmarkTimer("scanpy.read_h5ad [real pbmc 1k]", ITERATIONS)

        def _read():
            adata = sc.read_h5ad(REAL_H5AD)
            result = adata.shape
            del adata
            gc.collect()
            return result

        timer.run(_read)
        mem = MemoryProfiler.measure(lambda: sc.read_h5ad(REAL_H5AD))
        stats = timer.stats()

        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

        assert stats["mean_s"] < 120, "read_h5ad real should complete within 120s"


# ===================================================================
# organize_column_dtypes benchmarks
# ===================================================================

@pytest.mark.benchmark
class TestBenchOrganizeColumnDtypes:

    @staticmethod
    def _load_obs(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)
        obs = adata.obs.copy()
        del adata
        gc.collect()
        return obs

    def test_bench_organize_column_dtypes_small(self, small_h5ad):
        obs = self._load_obs(small_h5ad)
        timer = BenchmarkTimer("organize_column_dtypes [small]", ITERATIONS)
        timer.run(organize_column_dtypes, obs)
        mem = MemoryProfiler.measure(organize_column_dtypes, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    def test_bench_organize_column_dtypes_large(self, large_h5ad):
        obs = self._load_obs(large_h5ad)
        timer = BenchmarkTimer("organize_column_dtypes [large]", ITERATIONS)
        timer.run(organize_column_dtypes, obs)
        mem = MemoryProfiler.measure(organize_column_dtypes, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_organize_column_dtypes_real(self):
        obs = self._load_obs(REAL_H5AD)
        timer = BenchmarkTimer("organize_column_dtypes [real]", ITERATIONS)
        timer.run(organize_column_dtypes, obs)
        mem = MemoryProfiler.measure(organize_column_dtypes, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")


# ===================================================================
# get_annotation_columns benchmarks
# ===================================================================

@pytest.mark.benchmark
class TestBenchGetAnnotationColumns:

    @staticmethod
    def _load_obs(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)
        obs = adata.obs.copy()
        del adata
        gc.collect()
        return obs

    def test_bench_get_annotation_columns_small(self, small_h5ad):
        obs = self._load_obs(small_h5ad)
        timer = BenchmarkTimer("get_annotation_columns [small]", ITERATIONS)
        timer.run(get_annotation_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))

    def test_bench_get_annotation_columns_large(self, large_h5ad):
        obs = self._load_obs(large_h5ad)
        timer = BenchmarkTimer("get_annotation_columns [large]", ITERATIONS)
        timer.run(get_annotation_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_get_annotation_columns_real(self):
        obs = self._load_obs(REAL_H5AD)
        timer = BenchmarkTimer("get_annotation_columns [real]", ITERATIONS)
        timer.run(get_annotation_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))


# ===================================================================
# get_pseudotime_columns benchmarks
# ===================================================================

@pytest.mark.benchmark
class TestBenchGetPseudotimeColumns:

    @staticmethod
    def _load_obs(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)
        obs = adata.obs.copy()
        del adata
        gc.collect()
        return obs

    def test_bench_get_pseudotime_columns_small(self, small_h5ad):
        obs = self._load_obs(small_h5ad)
        timer = BenchmarkTimer("get_pseudotime_columns [small]", ITERATIONS)
        timer.run(get_pseudotime_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))

    def test_bench_get_pseudotime_columns_large(self, large_h5ad):
        obs = self._load_obs(large_h5ad)
        timer = BenchmarkTimer("get_pseudotime_columns [large]", ITERATIONS)
        timer.run(get_pseudotime_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_get_pseudotime_columns_real(self):
        obs = self._load_obs(REAL_H5AD)
        timer = BenchmarkTimer("get_pseudotime_columns [real]", ITERATIONS)
        timer.run(get_pseudotime_columns, obs)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats))


# ===================================================================
# Full pipeline simulation (mirrors /columns and /clusters endpoints)
# ===================================================================

@pytest.mark.benchmark
class TestBenchFullPipeline:
    """
    Simulates the exact code path of the /columns and /clusters endpoints
    without HTTP/DB overhead — pure function-level measurement.
    """

    @staticmethod
    def _columns_pipeline(filepath):
        """Mirrors files.py h5ad_columns endpoint logic."""
        adata = sc.read_h5ad(filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        anno = get_annotation_columns(adata.obs)
        pseudo = get_pseudotime_columns(adata.obs)
        result = {"anno_columns": anno, "pseudo_columns": pseudo}
        del adata
        gc.collect()
        return result

    @staticmethod
    def _clusters_pipeline(filepath, anno_column):
        """Mirrors files.py h5ad_cluster endpoint logic."""
        adata = sc.read_h5ad(filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        clusters = list(map(str, adata.obs[anno_column].value_counts().index))
        result = {"clusters": clusters}
        del adata
        gc.collect()
        return result

    # --- /columns pipeline ---

    def test_bench_columns_pipeline_small(self, small_h5ad):
        timer = BenchmarkTimer("columns_pipeline [small]", ITERATIONS)
        timer.run(self._columns_pipeline, small_h5ad)
        mem = MemoryProfiler.measure(self._columns_pipeline, small_h5ad)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    def test_bench_columns_pipeline_large(self, large_h5ad):
        timer = BenchmarkTimer("columns_pipeline [large]", ITERATIONS)
        timer.run(self._columns_pipeline, large_h5ad)
        mem = MemoryProfiler.measure(self._columns_pipeline, large_h5ad)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_columns_pipeline_real(self):
        timer = BenchmarkTimer("columns_pipeline [real pbmc 1k]", ITERATIONS)
        timer.run(self._columns_pipeline, REAL_H5AD)
        mem = MemoryProfiler.measure(self._columns_pipeline, REAL_H5AD)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    # --- /clusters pipeline ---

    def test_bench_clusters_pipeline_small(self, small_h5ad):
        timer = BenchmarkTimer("clusters_pipeline [small]", ITERATIONS)
        timer.run(self._clusters_pipeline, small_h5ad, "cluster")
        mem = MemoryProfiler.measure(self._clusters_pipeline, small_h5ad, "cluster")
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    def test_bench_clusters_pipeline_large(self, large_h5ad):
        timer = BenchmarkTimer("clusters_pipeline [large]", ITERATIONS)
        timer.run(self._clusters_pipeline, large_h5ad, "cluster")
        mem = MemoryProfiler.measure(self._clusters_pipeline, large_h5ad, "cluster")
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")

    @pytest.mark.skipif(REAL_H5AD is None, reason="pbmc_light_1000.h5ad not found")
    def test_bench_clusters_pipeline_real(self):
        """Use first annotation column found in real data."""
        adata = sc.read_h5ad(REAL_H5AD)
        obs = organize_column_dtypes(adata.obs)
        anno_cols = get_annotation_columns(obs)
        del adata
        gc.collect()

        anno_col = anno_cols[0] if anno_cols else "cluster"

        timer = BenchmarkTimer("clusters_pipeline [real pbmc 1k]", ITERATIONS)
        timer.run(self._clusters_pipeline, REAL_H5AD, anno_col)
        mem = MemoryProfiler.measure(self._clusters_pipeline, REAL_H5AD, anno_col)
        stats = timer.stats()
        print("\n" + format_benchmark_result(stats, mem))
        print(f"  JSON: {json.dumps({**stats, **{k: v for k, v in mem.items() if k != 'result'}})}")
