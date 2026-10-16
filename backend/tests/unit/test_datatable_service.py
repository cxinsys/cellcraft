"""
Unit tests for the datatable service layer (app/datatable/service.py) extracted
in PR-8.

Scope: pin the load/cache/format orchestration with the path resolver, cache,
and ``RemoteDataTable`` engine mocked. The terminal catch-all raises
``ValidationFailedError`` (-> 400) with ``str(e)`` as the detail — matching the
router's previous ``HTTPException(400, str(e))`` behaviour exactly.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.core.exceptions import ValidationFailedError
from app.datatable import service


def _vgt(file_name="d.h5ad", source="local"):
    v = MagicMock()
    v.file_name = file_name
    v.source = source
    return v


def _user(username="dtsvc"):
    u = MagicMock()
    u.username = username
    return u


class TestLoadData:
    def test_resolver_failure_wraps_to_400(self):
        with patch("app.file.path_resolver.resolve_data_file_path",
                   side_effect=RuntimeError("bad source")):
            with pytest.raises(ValidationFailedError) as exc:
                service.load_data(vgt_info=_vgt(), current_user=_user())
        assert exc.value.status_code == 400
        assert exc.value.detail == "bad source"

    def test_success_returns_vgt_payload(self):
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/p/f"), \
                patch("app.datatable.cache.get_cached_dataframe", return_value=MagicMock()), \
                patch("app.datatable.service.RemoteDataTable") as RT, \
                patch("app.datatable.service.transform_df_to_vgt_format",
                      return_value={"columns": ["c"], "rows": [{"c": 1}]}):
            RT.return_value.get_filtered_sorted_paginated_data.return_value = {
                "df": MagicMock(), "totalRecords": 5
            }
            out = service.load_data(vgt_info=_vgt(), current_user=_user())
        assert out == {
            "status": "success",
            "columns": ["c"],
            "rows": [{"c": 1}],
            "totalRecords": 5,
        }

    def test_cache_miss_loads_and_caches(self):
        with patch("app.file.path_resolver.resolve_data_file_path", return_value="/p/f"), \
                patch("app.datatable.cache.get_cached_dataframe", return_value=None), \
                patch("app.datatable.cache.set_cached_dataframe") as set_cache, \
                patch("app.datatable.service.load_tab_file", return_value=MagicMock()), \
                patch("app.datatable.service.RemoteDataTable") as RT, \
                patch("app.datatable.service.transform_df_to_vgt_format",
                      return_value={"columns": [], "rows": []}):
            RT.return_value.get_filtered_sorted_paginated_data.return_value = {
                "df": MagicMock(), "totalRecords": 0
            }
            service.load_data(vgt_info=_vgt(), current_user=_user())
        set_cache.assert_called_once()
