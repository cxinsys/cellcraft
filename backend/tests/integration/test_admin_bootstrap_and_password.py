"""
Initial admin provisioning, password change and forced change flow (#138).
"""
import os
import stat

import pytest

from app import models
from app.cli import create_admin, reset_password
from app.core import bootstrap
from app.core.config import settings
from app.core.exceptions import ValidationFailedError
from app.core.security import get_password_hash, verify_password

LEGACY_PASSWORD = bootstrap.LEGACY_DEFAULT_ADMIN_PASSWORDS[0]


@pytest.fixture
def password_file(tmp_path, monkeypatch):
    path = tmp_path / ".initial_admin_password"
    monkeypatch.setattr(settings, "INITIAL_ADMIN_PASSWORD_FILE", str(path))
    monkeypatch.setattr(settings, "CELLCRAFT_ADMIN_PASSWORD", None)
    monkeypatch.setattr(settings, "CELLCRAFT_ADMIN_EMAIL", "cellcraft@cellcraft.com")
    return path


def _user(db, email, password, superuser=False, username="u", must_change=False):
    user = models.User(username=username, email=email, hashed_password=get_password_hash(password),
                       is_active=True, is_superuser=superuser, must_change_password=must_change)
    db.add(user)
    db.commit()
    db.refresh(user)
    return user


def _login(client, email, password):
    return client.post("/routes/auth/login/access-token", data={"username": email, "password": password})


def _headers(client, email, password):
    return {"Authorization": f"Bearer {_login(client, email, password).json()['access_token']}"}


# ---------------------------------------------------------------------------
# bootstrap
# ---------------------------------------------------------------------------

class TestBootstrapAdmin:
    def test_new_install_creates_admin_with_generated_password(self, db_session, password_file, capsys):
        admin = bootstrap.bootstrap_admin(db_session)

        assert admin.is_superuser and admin.is_active and admin.must_change_password
        assert admin.username == "admin" and admin.email == "cellcraft@cellcraft.com"
        content = password_file.read_text()
        password = content.split("password: ")[1].strip()
        assert verify_password(password, admin.hashed_password)
        assert not verify_password(LEGACY_PASSWORD, admin.hashed_password)
        assert stat.S_IMODE(os.stat(password_file).st_mode) == 0o600
        # 비밀번호는 파일로만 안내하고 기동 로그에는 남기지 않는다
        output = capsys.readouterr().out
        assert password not in output and str(password_file) in output

    def test_password_printed_only_when_file_cannot_be_written(self, db_session, password_file, monkeypatch, capsys):
        def fail(*args, **kwargs):
            raise OSError("read-only file system")
        monkeypatch.setattr(bootstrap, "write_initial_password_file", fail)

        admin = bootstrap.bootstrap_admin(db_session)

        output = capsys.readouterr().out
        password = output.split("initial admin password: ")[1].split()[0]
        assert verify_password(password, admin.hashed_password)
        assert admin.must_change_password

    def test_password_file_removed_when_commit_fails(self, db_session, password_file, monkeypatch):
        monkeypatch.setattr(db_session, "commit", lambda: (_ for _ in ()).throw(RuntimeError("db down")))
        with pytest.raises(RuntimeError):
            bootstrap.bootstrap_admin(db_session)
        assert not password_file.exists()

    def test_recreated_after_all_admins_deleted(self, db_session, password_file):
        first = bootstrap.bootstrap_admin(db_session)
        db_session.delete(first)
        db_session.commit()

        second = bootstrap.bootstrap_admin(db_session)
        assert second is not None and second.is_superuser

    @pytest.mark.parametrize("setting, value", [
        ("CELLCRAFT_ADMIN_PASSWORD", "short"),
        ("CELLCRAFT_ADMIN_EMAIL", "not-an-email"),
    ])
    def test_invalid_settings_reported_without_creating(self, db_session, password_file, monkeypatch, capsys,
                                                         setting, value):
        monkeypatch.setattr(settings, setting, value)
        assert bootstrap.bootstrap_admin(db_session) is None
        assert db_session.query(models.User).count() == 0
        assert "Initial admin not created" in capsys.readouterr().out

    def test_empty_env_values_use_defaults(self, db_session, password_file, monkeypatch):
        # compose가 미설정 변수를 빈 문자열로 전달하는 경우
        monkeypatch.setattr(settings, "CELLCRAFT_ADMIN_EMAIL", "")
        monkeypatch.setattr(settings, "CELLCRAFT_ADMIN_PASSWORD", "")
        admin = bootstrap.bootstrap_admin(db_session)
        assert admin.email == "cellcraft@cellcraft.com" and admin.must_change_password

    def test_stale_password_file_removed(self, db_session, password_file):
        _user(db_session, "cellcraft@cellcraft.com", "changed-password", superuser=True, username="admin")
        bootstrap.write_initial_password_file("cellcraft@cellcraft.com", "old-password")

        bootstrap.bootstrap_admin(db_session)
        assert not password_file.exists()

    def test_password_file_kept_until_changed(self, db_session, password_file):
        bootstrap.bootstrap_admin(db_session)
        bootstrap.bootstrap_admin(db_session)  # 재기동: 아직 변경 전이면 유지
        assert password_file.exists()

    def test_password_from_env_not_forced_to_change(self, db_session, password_file, monkeypatch):
        monkeypatch.setattr(settings, "CELLCRAFT_ADMIN_PASSWORD", "operator-chosen-pw")
        admin = bootstrap.bootstrap_admin(db_session)

        assert verify_password("operator-chosen-pw", admin.hashed_password)
        assert not admin.must_change_password
        assert not password_file.exists()

    def test_existing_admin_kept_and_legacy_password_flagged(self, db_session, password_file):
        legacy = _user(db_session, "cellcraft@cellcraft.com", LEGACY_PASSWORD, superuser=True, username="admin")
        other = _user(db_session, "other-admin@example.com", "a-strong-password", superuser=True)

        assert bootstrap.bootstrap_admin(db_session) is None
        db_session.refresh(legacy)
        db_session.refresh(other)
        assert legacy.must_change_password and verify_password(LEGACY_PASSWORD, legacy.hashed_password)
        assert not other.must_change_password
        assert db_session.query(models.User).filter_by(is_superuser=True).count() == 2
        assert not password_file.exists()

    def test_not_recreated_when_admin_username_missing(self, db_session, password_file):
        _user(db_session, "renamed@example.com", "a-strong-password", superuser=True, username="root")
        assert bootstrap.bootstrap_admin(db_session) is None
        assert db_session.query(models.User).filter_by(username="admin").count() == 0

    def test_email_taken_by_regular_user(self, db_session, password_file):
        _user(db_session, "cellcraft@cellcraft.com", "a-strong-password")
        assert bootstrap.bootstrap_admin(db_session) is None
        assert db_session.query(models.User).filter_by(is_superuser=True).count() == 0


class TestPasswordPolicy:
    @pytest.mark.parametrize("password", ["short", "", None, LEGACY_PASSWORD])
    def test_rejected(self, password):
        with pytest.raises(ValidationFailedError):
            bootstrap.validate_new_password(password)

    def test_same_as_current_rejected(self):
        with pytest.raises(ValidationFailedError):
            bootstrap.validate_new_password("same-password", get_password_hash("same-password"))

    def test_accepted(self):
        bootstrap.validate_new_password("a-new-password", get_password_hash("old-password"))


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------

class TestPasswordChangeApi:
    def test_change_password(self, client, db_session, sample_user, auth_headers):
        response = client.post("/routes/auth/password", headers=auth_headers,
                               json={"current_password": "testpassword123", "new_password": "new-password-1"})

        assert response.status_code == 200
        assert response.json()["must_change_password"] is False
        assert _login(client, sample_user.email, "new-password-1").status_code == 200
        assert _login(client, sample_user.email, "testpassword123").status_code == 400

    @pytest.mark.parametrize("payload", [
        {"current_password": "wrong-password", "new_password": "new-password-1"},
        {"current_password": "testpassword123", "new_password": "short"},
        {"current_password": "testpassword123", "new_password": "testpassword123"},
    ])
    def test_rejected(self, client, auth_headers, payload):
        assert client.post("/routes/auth/password", headers=auth_headers, json=payload).status_code == 400

    def test_requires_auth(self, client):
        response = client.post("/routes/auth/password", json={"current_password": "a", "new_password": "b" * 8})
        assert response.status_code == 401

    def test_plugins_endpoint_no_longer_changes_password(self, client, sample_user, auth_headers):
        response = client.post("/routes/auth/plugins", headers=auth_headers, json={"password": "new-password-1"})
        assert response.status_code == 400
        assert _login(client, sample_user.email, "testpassword123").status_code == 200


class TestForcedPasswordChange:
    @pytest.fixture
    def flagged_user(self, db_session, password_file):
        user = _user(db_session, "cellcraft@cellcraft.com", "temp-password-1", superuser=True,
                     username="admin", must_change=True)
        bootstrap.write_initial_password_file(user.email, "temp-password-1")
        return user

    def test_login_reports_flag(self, client, flagged_user):
        response = _login(client, flagged_user.email, "temp-password-1")
        assert response.status_code == 200
        assert response.json()["user_info"]["must_change_password"] is True

    def test_other_api_blocked_until_changed(self, client, flagged_user, password_file):
        headers = _headers(client, flagged_user.email, "temp-password-1")

        blocked = client.get("/routes/workflow/me", headers=headers)
        assert blocked.status_code == 403
        assert blocked.json()["detail"]["code"] == "password_change_required"
        me = client.get("/routes/auth/me", headers=headers)
        assert me.status_code == 200 and me.json()["must_change_password"] is True

        changed = client.post("/routes/auth/password", headers=headers,
                              json={"current_password": "temp-password-1", "new_password": "new-password-1"})
        assert changed.status_code == 200
        assert client.get("/routes/workflow/me", headers=headers).status_code != 403
        assert not password_file.exists()

    def test_regular_user_not_affected(self, client, auth_headers):
        assert client.get("/routes/auth/me", headers=auth_headers).json()["must_change_password"] is False


# ---------------------------------------------------------------------------
# operator CLI
# ---------------------------------------------------------------------------

class TestCli:
    def test_create_admin_generated_password(self, db_session):
        password = create_admin(db_session, email="ops@example.com", username="ops")
        admin = db_session.query(models.User).filter_by(email="ops@example.com").one()
        assert admin.is_superuser and admin.must_change_password
        assert verify_password(password, admin.hashed_password)

    def test_create_admin_duplicate_email(self, db_session, sample_user):
        with pytest.raises(ValidationFailedError):
            create_admin(db_session, email=sample_user.email, username="x")

    def test_reset_password(self, db_session, sample_user):
        password = reset_password(db_session, email=sample_user.email)
        db_session.refresh(sample_user)
        assert sample_user.must_change_password
        assert verify_password(password, sample_user.hashed_password)

    def test_reset_password_keeps_inactive_unless_requested(self, db_session, sample_inactive_user):
        reset_password(db_session, email=sample_inactive_user.email)
        db_session.refresh(sample_inactive_user)
        assert not sample_inactive_user.is_active

        reset_password(db_session, email=sample_inactive_user.email, activate=True)
        db_session.refresh(sample_inactive_user)
        assert sample_inactive_user.is_active

    def test_reset_password_removes_initial_password_file(self, db_session, password_file):
        admin = bootstrap.bootstrap_admin(db_session)
        reset_password(db_session, email=admin.email)
        assert not password_file.exists()

    def test_create_admin_invalid_email(self, db_session):
        with pytest.raises(ValidationFailedError):
            create_admin(db_session, email="not-an-email", username="x")

    def test_reset_password_unknown_user(self, db_session):
        with pytest.raises(ValidationFailedError):
            reset_password(db_session, email="nobody@example.com")
