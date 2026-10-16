"""
Auth domain service layer.

Extracted from ``auth/router.py`` in PR-8 (Phase 3d). Holds the login / token
issuance business logic; the router keeps only the OAuth2 form dependency and
the ``response_model`` wiring.

Exception policy (PR-8): credential failures raise ``ValidationFailedError``
(-> 400), which the global handler in ``app.main`` maps onto the exact
``{"detail": ...}`` wire format (status code + detail string unchanged).
"""
from datetime import timedelta
from typing import Any

from sqlalchemy.orm import Session

from app.core.bootstrap import remove_initial_password_file, validate_new_password
from app.core.exceptions import ValidationFailedError
from app.core.security import create_access_token, get_password_hash, verify_password
from app.core.config import settings
from app.user import crud as crud_user
from app.user import service as user_service


def login(*, db: Session, form_data) -> Any:
    """Authenticate an OAuth2 form login and issue a JWT + user_info payload."""
    user = crud_user.authenticate(
        db, email=form_data.username, password=form_data.password
    )
    if not user:
        raise ValidationFailedError("Incorrect email or password")
    elif not crud_user.is_active(user):
        raise ValidationFailedError("please login to activate")
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    return {
        "access_token": create_access_token(
            user.id, expires_delta=access_token_expires
        ),
        "token_type": "bearer",
        "user_info": {
            "is_superuser": crud_user.is_superuser(user),
            "email": user.email,
            "is_active": user.is_active,
            "must_change_password": bool(user.must_change_password),
        }
    }


def change_password(*, db: Session, current_user, current_password: str, new_password: str) -> Any:
    """현재 비밀번호를 확인한 뒤 비밀번호를 변경하고 변경 필요 상태를 해제한다."""
    if not verify_password(current_password, current_user.hashed_password):
        raise ValidationFailedError("Incorrect current password")
    validate_new_password(new_password, current_user.hashed_password)

    current_user.hashed_password = get_password_hash(new_password)
    current_user.must_change_password = False
    db.add(current_user)
    db.commit()
    db.refresh(current_user)
    remove_initial_password_file(current_user.email)
    return user_service.read_me(current_user=current_user)
