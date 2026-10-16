"""
Unit tests for the user service layer (app/user/service.py) extracted in PR-8.

Scope: pin the registration / profile-read / update logic with the user crud +
filesystem boundary mocked. A duplicate-email signup raises
``ValidationFailedError`` (-> 400) with the frozen detail string.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.core.exceptions import ValidationFailedError
from app.user import service


def _user_in(email="a@b.com", username="alice"):
    u = MagicMock()
    u.email = email
    u.username = username
    return u


class TestRegister:
    def test_duplicate_email_raises_400(self):
        with patch("app.user.service.crud_user.get_user_by_email", return_value=MagicMock()):
            with pytest.raises(ValidationFailedError) as exc:
                service.register(db=MagicMock(), user_in=_user_in())
        assert exc.value.status_code == 400
        assert exc.value.detail == "Email already registered"

    def test_success_creates_user_and_folder(self):
        created = MagicMock()
        with patch("app.user.service.crud_user.get_user_by_email", return_value=None), \
                patch("app.user.service.crud_user.create_user", return_value=created), \
                patch("app.user.service.os.makedirs") as mk:
            out = service.register(db=MagicMock(), user_in=_user_in(username="bob"))
        assert out is created
        mk.assert_called_once_with("./user/bob/data", exist_ok=True)


class TestProfile:
    def test_read_me_shape(self):
        cu = MagicMock()
        cu.id, cu.email, cu.username, cu.is_superuser = 1, "a@b.com", "alice", False
        out = service.read_me(current_user=cu)
        assert out == {"id": 1, "email": "a@b.com", "username": "alice", "is_superuser": False}

    def test_update_plugins_delegates_to_crud(self):
        cu = MagicMock()
        cu.id = 5
        user_in = MagicMock()
        with patch("app.user.service.crud_user.update_user", return_value="updated") as upd:
            out = service.update_plugins(db=MagicMock(), current_user=cu, user_in=user_in)
        assert out == "updated"
        assert upd.call_args[1]["user_id"] == 5
