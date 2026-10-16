"""Initial admin account provisioning and password policy (#138).

Replaces the admin user that the plugin CSV sync used to create with a fixed,
publicly known password:

- new installs: an admin is created only when no superuser exists, with
  ``CELLCRAFT_ADMIN_PASSWORD`` or a random password announced once via
  ``INITIAL_ADMIN_PASSWORD_FILE`` (mode 0600) and the startup log; a generated
  password must be changed at first login
- existing installs: accounts are kept as they are; superusers still using the
  legacy default password are flagged ``must_change_password``
"""
import logging
import os
import secrets
from typing import Optional

from sqlalchemy.orm import Session

from app.core.config import settings
from app.core.exceptions import ValidationFailedError
from app.core.security import get_password_hash, verify_password
from app.user import models

logger = logging.getLogger(__name__)

DEFAULT_ADMIN_USERNAME = "admin"
MIN_PASSWORD_LENGTH = 8
# 이전 버전이 모든 설치에 동일하게 생성하던 관리자 비밀번호 (기존 설치 감지·재사용 금지 용도)
LEGACY_DEFAULT_ADMIN_PASSWORDS = ("cellcraft2024!",)


def generate_password() -> str:
    return secrets.token_urlsafe(12)


def is_legacy_default_password(password: str) -> bool:
    return password in LEGACY_DEFAULT_ADMIN_PASSWORDS


def uses_legacy_default_password(user: models.User) -> bool:
    return any(verify_password(p, user.hashed_password) for p in LEGACY_DEFAULT_ADMIN_PASSWORDS)


def validate_new_password(new_password: str, current_hashed_password: Optional[str] = None) -> None:
    """비밀번호 규칙: 8자 이상, 현재 비밀번호와 다름, 이전 기본 비밀번호 금지."""
    if not isinstance(new_password, str) or len(new_password) < MIN_PASSWORD_LENGTH:
        raise ValidationFailedError(f"Password must be at least {MIN_PASSWORD_LENGTH} characters")
    if is_legacy_default_password(new_password):
        raise ValidationFailedError("This password cannot be used")
    if current_hashed_password and verify_password(new_password, current_hashed_password):
        raise ValidationFailedError("New password must be different from the current password")


def write_initial_password_file(email: str, password: str, path: Optional[str] = None) -> str:
    path = path or settings.INITIAL_ADMIN_PASSWORD_FILE
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as f:
        f.write(f"email: {email}\npassword: {password}\n")
    os.chmod(path, 0o600)
    return path


def remove_initial_password_file(email: str, path: Optional[str] = None) -> bool:
    """초기 비밀번호 파일이 해당 계정의 것이면 삭제한다."""
    path = path or settings.INITIAL_ADMIN_PASSWORD_FILE
    try:
        with open(path) as f:
            if f"email: {email}\n" not in f.read():
                return False
        os.remove(path)
        return True
    except OSError:
        return False


def bootstrap_admin(db: Session) -> Optional[models.User]:
    """superuser가 없으면 초기 관리자를 만들고, 기본 비밀번호를 쓰는 관리자를 표시한다.

    Returns:
        새로 만든 관리자 (없으면 None)
    """
    superusers = db.query(models.User).filter(models.User.is_superuser.is_(True)).all()
    created = None

    if not superusers:
        email = settings.CELLCRAFT_ADMIN_EMAIL
        if db.query(models.User).filter_by(email=email).first():
            logger.warning(f"Initial admin not created: email {email} is already used by a non-admin account")
            return None
        password = settings.CELLCRAFT_ADMIN_PASSWORD
        generated = not password
        if generated:
            password = generate_password()
        else:
            validate_new_password(password)
        created = models.User(
            username=DEFAULT_ADMIN_USERNAME,
            email=email,
            hashed_password=get_password_hash(password),
            is_active=True,
            is_superuser=True,
            must_change_password=generated,
        )
        db.add(created)
        db.commit()
        if generated:
            path = write_initial_password_file(email, password)
            banner = "=" * 80
            print(f"{banner}\nInitial admin account created\n  email:    {email}\n  password: {password}\n"
                  f"  (saved to {os.path.abspath(path)}; change it at first login)\n{banner}", flush=True)
        else:
            print(f"Initial admin account created: {email} (password from CELLCRAFT_ADMIN_PASSWORD)", flush=True)
        return created

    flagged = 0
    for user in superusers:
        if not user.must_change_password and uses_legacy_default_password(user):
            user.must_change_password = True
            flagged += 1
    if flagged:
        db.commit()
        print(f"WARNING: {flagged} admin account(s) use the legacy default password; "
              "a password change is required at next login", flush=True)
    return created
