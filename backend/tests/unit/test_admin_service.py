"""
Unit tests for the admin service layer (app/admin/service.py) extracted in PR-8.

Scope: pin the admin-gate + listing/mutation business logic with the admin crud,
Docker, and plugin-sync boundaries mocked. Asserts the domain-exception mapping:
the superuser gate raises ``PermissionDeniedError`` (-> 403) and the "not found"
guards raise ``NotFoundError`` (-> 404), which the global handler renders as the
unchanged ``{"detail": ...}`` body.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.core.exceptions import (
    NotFoundError, PermissionDeniedError, ExternalServiceError
)
from app.admin import service


def _admin():
    u = MagicMock()
    u.is_superuser = True
    return u


def _nonadmin():
    u = MagicMock()
    u.is_superuser = False
    return u


_LIST_KW = dict(amount=10, page_num=1, sort="id", order="asc", searchTerm="")


# ---------------------------------------------------------------------------
# admin gate
# ---------------------------------------------------------------------------

class TestAdminGate:
    def test_nonadmin_blocked_403(self):
        with pytest.raises(PermissionDeniedError) as exc:
            service.get_users_count(db=MagicMock(), current_user=_nonadmin())
        assert exc.value.status_code == 403
        assert exc.value.detail == "Access denied: Admins only"

    def test_nonadmin_blocked_on_listing(self):
        with pytest.raises(PermissionDeniedError):
            service.get_filtered_users(db=MagicMock(), current_user=_nonadmin(), **_LIST_KW)


# ---------------------------------------------------------------------------
# listings + counts
# ---------------------------------------------------------------------------

class TestListings:
    def test_users_empty_raises_404(self):
        with patch("app.admin.service.crud_admin.get_filtered_users", return_value=([], 0)):
            with pytest.raises(NotFoundError) as exc:
                service.get_filtered_users(db=MagicMock(), current_user=_admin(), **_LIST_KW)
        assert exc.value.status_code == 404
        assert exc.value.detail == "Users not found"

    def test_users_returns_rows(self):
        rows = [MagicMock(), MagicMock()]
        with patch("app.admin.service.crud_admin.get_filtered_users", return_value=(rows, 2)):
            out = service.get_filtered_users(db=MagicMock(), current_user=_admin(), **_LIST_KW)
        assert out is rows

    def test_users_count_returns_value(self):
        with patch("app.admin.service.crud_admin.get_users_count", return_value=42):
            assert service.get_users_count(db=MagicMock(), current_user=_admin()) == 42

    def test_files_formats_rows(self):
        f = MagicMock()
        f.id, f.file_name, f.file_path, f.file_size, f.folder, f.user_id = 1, "a.h5ad", "/p", 10, "fold", 5
        with patch("app.admin.service.crud_admin.get_filtered_files", return_value=([(f, "bob")], 1)):
            out = service.get_filtered_files(db=MagicMock(), current_user=_admin(), **_LIST_KW)
        assert out["total_count"] == 1
        assert out["data"][0]["username"] == "bob"
        assert out["data"][0]["file_name"] == "a.h5ad"

    def test_files_empty_raises_404(self):
        with patch("app.admin.service.crud_admin.get_filtered_files", return_value=([], 0)):
            with pytest.raises(NotFoundError) as exc:
                service.get_filtered_files(db=MagicMock(), current_user=_admin(), **_LIST_KW)
        assert exc.value.detail == "Files not found"


# ---------------------------------------------------------------------------
# mutations
# ---------------------------------------------------------------------------

class TestMutations:
    def test_delete_user_missing_raises_404(self):
        with patch("app.admin.service.crud_admin.delete_user", return_value=False):
            with pytest.raises(NotFoundError) as exc:
                service.delete_user(db=MagicMock(), current_user=_admin(), user_id=1)
        assert exc.value.detail == "User not found"

    def test_delete_user_success_message(self):
        with patch("app.admin.service.crud_admin.delete_user", return_value=True):
            out = service.delete_user(db=MagicMock(), current_user=_admin(), user_id=1)
        assert out == {"message": "User deleted successfully"}

    def test_update_user_missing_raises_404(self):
        with patch("app.admin.service.crud_admin.update_user", return_value=None):
            with pytest.raises(NotFoundError):
                service.update_user(db=MagicMock(), current_user=_admin(), user_id=1, user_data={})

    def test_cancel_task_missing_raises_404(self):
        with patch("app.admin.service.crud_admin.cancel_task", return_value=False):
            with pytest.raises(NotFoundError) as exc:
                service.cancel_task(db=MagicMock(), current_user=_admin(), task_id=9)
        assert exc.value.detail == "Task not found"


# ---------------------------------------------------------------------------
# plugin sync / container manager
# ---------------------------------------------------------------------------

class TestPluginSync:
    def test_sync_status_success(self):
        mgr = MagicMock()
        mgr.get_sync_status.return_value = {"ok": True}
        with patch("app.admin.service.PluginSyncManager", return_value=mgr):
            out = service.get_plugin_sync_status(current_user=_admin())
        assert out == {"ok": True}

    def test_sync_status_failure_raises_external(self):
        with patch("app.admin.service.PluginSyncManager", side_effect=RuntimeError("git down")):
            with pytest.raises(ExternalServiceError) as exc:
                service.get_plugin_sync_status(current_user=_admin())
        assert exc.value.status_code == 500
        assert "Failed to get sync status" in exc.value.detail

    def test_container_status_success(self):
        with patch("app.admin.service.container_manager") as cm:
            cm.get_status.return_value = {"tracked_tasks": 0}
            out = service.get_container_manager_status(current_user=_admin())
        assert out == {"tracked_tasks": 0}
