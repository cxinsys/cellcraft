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

from app.core.exceptions import ValidationFailedError
from app.core.security import create_access_token
from app.core.config import settings
from app.user import crud as crud_user


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
            "is_active": user.is_active
        }
    }
