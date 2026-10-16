"""
User domain service layer.

Extracted from ``auth/router.py`` in PR-8 (Phase 3d). Holds the user registration
/ profile-read / plugin-update business logic. The three endpoints physically
remain on ``auth/router.py`` (moving them would change their mounted paths and
break the OpenAPI contract snapshot), but their logic now lives here.

Exception policy (PR-8): a duplicate-email signup raises ``ValidationFailedError``
(-> 400); the global handler maps it onto the exact ``{"detail": ...}`` wire
format (status code + detail string unchanged).
"""
import os
from typing import Any

from sqlalchemy.orm import Session

from app import models
from app.core.exceptions import ValidationFailedError
from app.user import crud as crud_user
from app.user import schemas as user_schemas


def register(*, db: Session, user_in: user_schemas.UserCreate) -> Any:
    """Create a new user, provisioning their data folder (400 on duplicate email)."""
    existing = crud_user.get_user_by_email(db, email=user_in.email)
    if existing:
        raise ValidationFailedError("Email already registered")
    created = crud_user.create_user(db, user=user_in)
    USER_DIRECTORY_NAME = './user/' + user_in.username + '/data'
    os.makedirs(USER_DIRECTORY_NAME, exist_ok=True)
    # 회원가입 시 보내는 확인 이메일
    # if user_in.email:
    #     send_new_account_email(
    #         email_to=user_in.email, username=user_in.email, password=user_in.password
    #     )
    return created


def read_me(*, current_user: models.User) -> Any:
    """Return the current user's public profile fields."""
    return {
        "id": current_user.id,
        "email": current_user.email,
        "username": current_user.username,
        "is_superuser": current_user.is_superuser,
    }


def update_plugins(*, db: Session, current_user: models.User,
                   user_in: user_schemas.UserUpdate) -> Any:
    """Update the current user's record (plugin associations / password)."""
    return crud_user.update_user(db, user_id=current_user.id, user=user_in)
