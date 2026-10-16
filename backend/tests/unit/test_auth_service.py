"""
Unit tests for the auth service layer (app/auth/service.py) extracted in PR-8.

Scope: pin the login / token-issuance logic with the user crud + token boundary
mocked. Credential failures raise ``ValidationFailedError`` (-> 400); the detail
strings are frozen to match the auth characterization tests.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.core.exceptions import ValidationFailedError
from app.auth import service


def _form(username="a@b.com", password="pw"):
    f = MagicMock()
    f.username = username
    f.password = password
    return f


class TestLogin:
    def test_bad_credentials_raises_400(self):
        with patch("app.auth.service.crud_user.authenticate", return_value=None):
            with pytest.raises(ValidationFailedError) as exc:
                service.login(db=MagicMock(), form_data=_form())
        assert exc.value.status_code == 400
        assert exc.value.detail == "Incorrect email or password"

    def test_inactive_user_raises_400(self):
        user = MagicMock()
        with patch("app.auth.service.crud_user.authenticate", return_value=user), \
                patch("app.auth.service.crud_user.is_active", return_value=False):
            with pytest.raises(ValidationFailedError) as exc:
                service.login(db=MagicMock(), form_data=_form())
        assert exc.value.status_code == 400
        assert exc.value.detail == "please login to activate"

    def test_success_returns_token_payload(self):
        user = MagicMock()
        user.id = 7
        user.email = "a@b.com"
        user.is_active = True
        with patch("app.auth.service.crud_user.authenticate", return_value=user), \
                patch("app.auth.service.crud_user.is_active", return_value=True), \
                patch("app.auth.service.crud_user.is_superuser", return_value=False), \
                patch("app.auth.service.create_access_token", return_value="tok"):
            out = service.login(db=MagicMock(), form_data=_form())
        assert set(out.keys()) == {"access_token", "token_type", "user_info"}
        assert out["access_token"] == "tok"
        assert out["token_type"] == "bearer"
        assert out["user_info"] == {"is_superuser": False, "email": "a@b.com", "is_active": True}
