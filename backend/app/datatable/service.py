"""
Datatable domain service layer.

Extracted from ``datatable/router.py`` in PR-8 (Phase 3d). Holds the h5ad/tabular
data-load orchestration: path resolution, dataframe caching, and the
filter/sort/paginate + vue-good-table formatting. The ``RemoteDataTable`` engine
and the workflow/cache utilities are invoked, never modified (utils internals are
out of scope per PR-8).

Exception policy (PR-8): the terminal catch-all now raises
``ValidationFailedError`` (-> 400); the global handler maps it onto the exact
``{"detail": str(e)}`` wire format, so the response (400 + detail string) is
unchanged.
"""
from typing import Any

from app import models
from app.core.exceptions import ValidationFailedError
from app.workflow.utils import load_tab_file, transform_df_to_vgt_format
from app.datatable.utils import RemoteDataTable, DataTableEvent


def load_data(*, vgt_info: DataTableEvent, current_user: models.User) -> Any:
    """Load, cache, filter/sort/paginate and vgt-format a user's data file."""
    try:
        from app.file.path_resolver import resolve_data_file_path
        from app.datatable.cache import get_cached_dataframe, set_cached_dataframe

        file_path = resolve_data_file_path(vgt_info.file_name, current_user.username, vgt_info.source)

        df = get_cached_dataframe(file_path)
        if df is None:
            df = load_tab_file(file_path)
            set_cached_dataframe(file_path, df)

        remote_table = RemoteDataTable(df, vgt_info)
        data = remote_table.get_filtered_sorted_paginated_data()
        vgt_data = transform_df_to_vgt_format(data["df"])
        return {"status": "success", "columns": vgt_data['columns'], "rows": vgt_data["rows"], "totalRecords": data["totalRecords"]}
    except Exception as e:
        raise ValidationFailedError(str(e))
