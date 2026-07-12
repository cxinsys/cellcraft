"""
Characterization tests for the auth login endpoint.

Purpose: pin the CURRENT behavior of ``POST /routes/auth/login/access-token``
exactly as it is today, so refactoring PRs (service extraction, etc.) can prove
they did not change observable behavior. These tests describe reality, not the
ideal — if something looks odd, it is preserved intentionally.

Source of truth read while writing these tests:
- app/routes/endpoints/auth.py :: login_access_token (lines 47-69)
- app/routes/api.py mounts the auth router; app/main.py mounts api_router with
  prefix settings.ROUTES_STR ("/routes"), and auth router is mounted at
  "/auth" -> full prefix "/routes/auth".

Pinned behavior:
- Successful login returns 200 with keys: access_token, token_type ("bearer"),
  and a nested user_info dict {is_superuser, email, is_active}.
- Wrong password / unknown user -> 400 "Incorrect email or password".
- Inactive user -> 400 "please login to activate".
- OAuth2 form uses the 'username' field to carry the email.
"""
import pytest
from fastapi.testclient import TestClient

from app.database import models
from app.core.config import settings

LOGIN_URL = f"{settings.ROUTES_STR}/auth/login/access-token"


@pytest.mark.integration
@pytest.mark.characterization
@pytest.mark.auth
class TestCharacterizationAuthLogin:
    """Freeze the current token issuance response shape and error handling."""

    def test_login_success_response_shape(
        self,
        client: TestClient,
        sample_user: models.User,
    ):
        """Successful login returns the exact current token payload shape."""
        response = client.post(
            LOGIN_URL,
            data={"username": sample_user.email, "password": "testpassword123"},
        )

        assert response.status_code == 200
        data = response.json()

        # Top-level keys are exactly these three (current contract).
        assert set(data.keys()) == {"access_token", "token_type", "user_info"}
        assert data["token_type"] == "bearer"
        assert isinstance(data["access_token"], str) and data["access_token"]

        # user_info is a nested dict with these three keys, in this shape.
        user_info = data["user_info"]
        assert set(user_info.keys()) == {"is_superuser", "email", "is_active"}
        assert user_info["email"] == sample_user.email
        assert user_info["is_active"] is True
        assert user_info["is_superuser"] is False

    def test_login_incorrect_password_returns_400(
        self,
        client: TestClient,
        sample_user: models.User,
    ):
        """Wrong password yields 400 with the current fixed detail message."""
        response = client.post(
            LOGIN_URL,
            data={"username": sample_user.email, "password": "wrongpassword"},
        )

        assert response.status_code == 400
        assert response.json()["detail"] == "Incorrect email or password"

    def test_login_nonexistent_user_returns_400(
        self,
        client: TestClient,
    ):
        """Unknown user yields the same 400 message as a wrong password."""
        response = client.post(
            LOGIN_URL,
            data={"username": "nobody@example.com", "password": "whatever"},
        )

        assert response.status_code == 400
        assert response.json()["detail"] == "Incorrect email or password"

    def test_login_inactive_user_rejected(
        self,
        client: TestClient,
        sample_inactive_user: models.User,
    ):
        """Inactive user (correct credentials) is rejected with 400 'activate'."""
        response = client.post(
            LOGIN_URL,
            data={"username": sample_inactive_user.email, "password": "testpassword123"},
        )

        assert response.status_code == 400
        assert response.json()["detail"] == "please login to activate"

    def test_login_missing_form_fields_returns_422(
        self,
        client: TestClient,
    ):
        """Missing OAuth2 form fields fail FastAPI validation with 422."""
        # No form body at all -> OAuth2PasswordRequestForm required fields missing.
        response = client.post(LOGIN_URL)
        assert response.status_code == 422
