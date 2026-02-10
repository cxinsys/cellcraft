"""H5AD file compression utility.

Compresses uploaded H5AD files using HDF5 internal gzip compression.
sc.read_h5ad() transparently decompresses, so no downstream code changes needed.
"""

import gc
import logging
import os
import shutil
import tempfile
from typing import Tuple

import scanpy as sc

logger = logging.getLogger(__name__)

# anndata 0.8.0 reserves '_index' as a column name.
# Some H5AD files (e.g., from older pipelines) have '_index' in raw.var or obs,
# which causes write_h5ad() to raise ValueError.
_ANNDATA_RESERVED_COLUMNS = {"_index"}


def _fix_reserved_columns(adata) -> None:
    """Rename reserved column names in-place so write_h5ad() succeeds."""
    for attr_name in ("obs", "var"):
        df = getattr(adata, attr_name, None)
        if df is not None:
            conflicts = _ANNDATA_RESERVED_COLUMNS & set(df.columns)
            if conflicts:
                rename_map = {c: f"{c}_original" for c in conflicts}
                df.rename(columns=rename_map, inplace=True)
                logger.info("Renamed reserved columns in %s: %s", attr_name, rename_map)

    # raw.var is accessed via adata.raw.var (Raw object)
    if adata.raw is not None:
        raw_var = adata.raw.var
        conflicts = _ANNDATA_RESERVED_COLUMNS & set(raw_var.columns)
        if conflicts:
            rename_map = {c: f"{c}_original" for c in conflicts}
            raw_var_fixed = raw_var.rename(columns=rename_map)
            # Reconstruct raw with fixed var (Raw is immutable, so recreate)
            import anndata
            adata._raw = anndata.Raw(
                adata,
                X=adata.raw.X,
                var=raw_var_fixed,
                varm=adata.raw.varm if hasattr(adata.raw, "varm") else None,
            )
            logger.info("Renamed reserved columns in raw.var: %s", rename_map)


def compress_h5ad_file(
    file_path: str,
    min_size_bytes: int = 1 * 1024 * 1024,
    compression: str = "gzip",
) -> Tuple[bool, int]:
    """Compress an H5AD file in-place using HDF5 internal compression.

    Args:
        file_path: Absolute path to the .h5ad file.
        min_size_bytes: Skip compression if file is smaller than this (default 1MB).
        compression: HDF5 compression algorithm (default 'gzip').

    Returns:
        (compressed, final_size): Whether compression was applied, and the final file size in bytes.
    """
    if not file_path.lower().endswith(".h5ad"):
        logger.debug("Skipping non-H5AD file: %s", file_path)
        return False, os.path.getsize(file_path)

    original_size = os.path.getsize(file_path)

    if original_size < min_size_bytes:
        logger.debug(
            "Skipping compression for %s (size %d < min %d)",
            file_path, original_size, min_size_bytes,
        )
        return False, original_size

    tmp_fd = None
    tmp_path = None
    adata = None

    try:
        # Create temp file in same directory to ensure same filesystem (atomic move)
        dir_name = os.path.dirname(file_path)
        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".h5ad.tmp", dir=dir_name)
        os.close(tmp_fd)
        tmp_fd = None

        adata = sc.read_h5ad(file_path)
        _fix_reserved_columns(adata)
        adata.write_h5ad(tmp_path, compression=compression)

        compressed_size = os.path.getsize(tmp_path)

        if compressed_size < original_size:
            shutil.move(tmp_path, file_path)
            tmp_path = None  # Moved successfully, don't clean up
            savings_pct = (1 - compressed_size / original_size) * 100
            logger.info(
                "Compressed %s: %d -> %d bytes (%.1f%% reduction)",
                file_path, original_size, compressed_size, savings_pct,
            )
            return True, compressed_size
        else:
            logger.info(
                "Compression not beneficial for %s (%d >= %d), keeping original",
                file_path, compressed_size, original_size,
            )
            return False, original_size

    except Exception:
        logger.warning(
            "Failed to compress %s, keeping original file",
            file_path, exc_info=True,
        )
        return False, original_size

    finally:
        if adata is not None:
            del adata
        gc.collect()

        # Clean up temp file if it still exists
        if tmp_path is not None and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        if tmp_fd is not None:
            try:
                os.close(tmp_fd)
            except OSError:
                pass
