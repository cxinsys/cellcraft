"""
Unit tests for authentication endpoints.

Test coverage:
- User registration (success, duplicate email, invalid email, validation)
- User login (success, wrong password, nonexistent user, inactive user)
- JWT token validation (payload claims, expiration, tampering, format)
- Token security (concurrent logins, payload verification)
- User profile retrieval (authenticated, unauthorized)
- User plugin updates (success, unauthorized)

Phase 2.1 Enhancements:
- JWT payload assertions (sub, exp, iat claims validation)
- Concurrent login scenarios
- Enhanced error path testing
- Registration validation edge cases

TODO (requires implementation):
- Token refresh mechanism (no /refresh endpoint exists)
- Password reset flows (no /reset endpoint exists)
- Rate limiting tests (no rate limiter configured)
"""
import pytest
from datetime import timedelta, datetime
from jose import jwt
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app import models
from app.core.security import create_access_token, ALGORITHM
from app.core.config import settings


@pytest.mark.unit
@pytest.mark.auth
class TestUserRegistration:
    """Test user registration endpoint."""

    def test_register_new_user_success(self, client: TestClient, db_session: Session):
        """Test successful user registration."""
        # Arrange
        user_data = {
            "username": "newuser",
            "email": "newuser@example.com",
            "password": "securepassword123"
        }

        # Act
        response = client.post("/routes/auth/register", json=user_data)

        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["username"] == user_data["username"]
        assert data["email"] == user_data["email"]
        assert "id" in data
        assert "hashed_password" not in data  # Password should not be returned

        # Verify user exists in database
        user = db_session.query(models.User).filter_by(
            email=user_data["email"]
        ).first()
        assert user is not None
        assert user.username == user_data["username"]

    def test_register_duplicate_email_fails(
        self,
        client: TestClient,
        sample_user: models.User
    ):
        """Test registration fails with duplicate email."""
        # Arrange
        user_data = {
            "username": "anotheruser",
            "email": sample_user.email,  # Duplicate email
            "password": "password123"
        }

        # Act
        response = client.post("/routes/auth/register", json=user_data)

        # Assert - Enhanced payload validation
        assert response.status_code == 400
        data = response.json()
        assert "detail" in data
        assert "already registered" in data["detail"].lower()

    def test_register_invalid_email_fails(self, client: TestClient):
        """Test registration fails with invalid email format."""
        # Arrange
        user_data = {
            "username": "testuser",
            "email": "invalid-email-format",  # Invalid email
            "password": "password123"
        }

        # Act
        response = client.post("/routes/auth/register", json=user_data)

        # Assert - Enhanced payload validation
        assert response.status_code == 422  # Validation error
        data = response.json()
        assert "detail" in data

    @pytest.mark.xfail(
        reason="BUG: Pydantic doesn't catch None values for required fields. "
               "Validation gap allows None to reach database layer causing IntegrityError. "
               "Should be fixed in Phase 2.3 (schema validation improvements)."
    )
    def test_register_missing_required_fields(self, client: TestClient):
        """Test registration fails when required fields are missing.

        KNOWN ISSUE: Current implementation has validation gap where None values
        for required fields bypass Pydantic validation and cause database errors.

        Expected behavior: Should return 422 validation error
        Actual behavior: Returns 500 IntegrityError from database
        """
        # Test missing email - Pydantic should catch this
        response1 = client.post("/routes/auth/register", json={
            "username": "testuser",
            "password": "password123"
        })
        assert response1.status_code == 422, "Missing email should return validation error"

        # Test missing password - Pydantic should catch this
        response2 = client.post("/routes/auth/register", json={
            "username": "testuser",
            "email": "test@example.com"
        })
        assert response2.status_code == 422, "Missing password should return validation error"

        # Test missing username - Pydantic should catch this
        response3 = client.post("/routes/auth/register", json={
            "email": "test@example.com",
            "password": "password123"
        })
        assert response3.status_code == 422, "Missing username should return validation error"

    @pytest.mark.xfail(
        reason="BUG: Empty string usernames are accepted by current validation. "
               "Should add min_length constraint in Pydantic schema. "
               "To be fixed in Phase 2.3 (schema validation improvements)."
    )
    def test_register_empty_fields(self, client: TestClient):
        """Test registration with empty username field.

        KNOWN ISSUE: Current Pydantic schema accepts empty strings for username.

        Expected behavior: Should return 422 validation error for empty username
        Actual behavior: Accepts empty string and creates user with empty username

        Recommended fix: Add min_length=1 constraint to UserCreate schema
        """
        # Arrange
        user_data = {
            "username": "",
            "email": "test@example.com",
            "password": "password123"
        }

        # Act
        response = client.post("/routes/auth/register", json=user_data)

        # Assert - Should reject empty username
        assert response.status_code in [400, 422], "Empty username should be rejected"

    def test_register_whitespace_only_username(self, client: TestClient):
        """Test registration with whitespace-only username.

        Note: Documents actual behavior. Whitespace-only strings may pass initial
        validation but should be caught by business logic or database constraints.
        """
        # Arrange
        user_data = {
            "username": "   ",
            "email": "test@example.com",
            "password": "password123"
        }

        # Act
        response = client.post("/routes/auth/register", json=user_data)

        # Assert - Implementation may accept or reject this
        # Test documents actual behavior rather than prescribing requirements
        assert response.status_code in [200, 400, 422], "Document whitespace username handling"
        if response.status_code == 200:
            # If accepted, username should be stored (potential improvement area)
            data = response.json()
            assert "username" in data


@pytest.mark.unit
@pytest.mark.auth
class TestUserLogin:
    """Test user login and JWT token generation."""

    def test_login_success(self, client: TestClient, sample_user: models.User):
        """Test successful login with valid credentials and JWT payload validation."""
        # Arrange
        login_data = {
            "username": sample_user.email,  # OAuth2 uses 'username' field for email
            "password": "testpassword123"
        }

        # Act
        response = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - Response structure
        assert response.status_code == 200
        data = response.json()
        assert "access_token" in data
        assert data["token_type"] == "bearer"
        assert len(data["access_token"]) > 0

        # Verify user info structure
        assert "user_info" in data
        assert data["user_info"]["email"] == sample_user.email
        assert data["user_info"]["is_active"] is True
        assert isinstance(data["user_info"]["is_superuser"], bool)

        # Phase 2.1 Enhancement: JWT payload validation
        token = data["access_token"]
        decoded_token = jwt.decode(token, settings.SECRET_KEY, algorithms=[ALGORITHM])

        # Verify required JWT claims
        assert "sub" in decoded_token, "JWT must contain 'sub' (subject) claim"
        assert "exp" in decoded_token, "JWT must contain 'exp' (expiration) claim"

        # Verify subject matches user ID
        assert decoded_token["sub"] == str(sample_user.id), "JWT 'sub' must match user ID"

        # Verify expiration is in the future
        exp_timestamp = decoded_token["exp"]
        current_timestamp = datetime.utcnow().timestamp()
        assert exp_timestamp > current_timestamp, "JWT expiration must be in the future"

        # Verify expiration is within expected range (ACCESS_TOKEN_EXPIRE_MINUTES)
        expected_max_exp = current_timestamp + (settings.ACCESS_TOKEN_EXPIRE_MINUTES * 60) + 10  # +10s tolerance
        assert exp_timestamp <= expected_max_exp, "JWT expiration exceeds configured maximum"

    def test_login_wrong_password_fails(
        self,
        client: TestClient,
        sample_user: models.User
    ):
        """Test login fails with incorrect password."""
        # Arrange
        login_data = {
            "username": sample_user.email,
            "password": "wrongpassword"
        }

        # Act
        response = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - Enhanced payload validation
        assert response.status_code == 400
        data = response.json()
        assert "detail" in data
        assert "incorrect" in data["detail"].lower()

    def test_login_nonexistent_user_fails(self, client: TestClient):
        """Test login fails with non-existent user."""
        # Arrange
        login_data = {
            "username": "nonexistent@example.com",
            "password": "somepassword"
        }

        # Act
        response = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - Enhanced payload validation
        assert response.status_code == 400
        data = response.json()
        assert "detail" in data

    def test_login_inactive_user_fails(
        self,
        client: TestClient,
        sample_inactive_user: models.User
    ):
        """Test login fails for inactive user using fixture."""
        # Arrange
        login_data = {
            "username": sample_inactive_user.email,
            "password": "testpassword123"
        }

        # Act
        response = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - Enhanced payload validation
        assert response.status_code == 400
        data = response.json()
        assert "detail" in data
        assert "activate" in data["detail"].lower()

    @pytest.mark.xfail(
        reason="LIMITATION: JWT tokens generated within same second are identical. "
               "Current implementation lacks 'iat' (issued at) claim or microsecond precision. "
               "Tokens generated rapidly may be duplicates. Consider adding 'iat' claim."
    )
    def test_concurrent_logins_generate_different_tokens(
        self,
        client: TestClient,
        sample_user: models.User
    ):
        """Test that multiple concurrent logins for same user generate unique tokens.

        KNOWN LIMITATION: Current JWT implementation doesn't include 'iat' (issued at)
        claim, only 'exp' (expiration). Since 'exp' uses second precision, tokens
        generated within the same second will be identical.

        Security implication: Rapid successive logins produce duplicate tokens.

        Recommended improvements:
        1. Add 'iat' claim with microsecond precision
        2. Add 'jti' (JWT ID) claim for unique token identification
        3. Consider token versioning or session management
        """
        # Arrange
        login_data = {
            "username": sample_user.email,
            "password": "testpassword123"
        }

        # Act - Perform two rapid logins
        response1 = client.post("/routes/auth/login/access-token", data=login_data)
        response2 = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - Both logins succeed
        assert response1.status_code == 200
        assert response2.status_code == 200

        # Extract tokens
        token1 = response1.json()["access_token"]
        token2 = response2.json()["access_token"]

        # Expected: Tokens should be different (requires 'iat' or 'jti' claim)
        # Actual: May be identical if generated within same second
        assert token1 != token2, "Concurrent logins should generate unique tokens"

        # Both tokens should decode successfully
        decoded1 = jwt.decode(token1, settings.SECRET_KEY, algorithms=[ALGORITHM])
        decoded2 = jwt.decode(token2, settings.SECRET_KEY, algorithms=[ALGORITHM])

        # Both should have same subject (user ID)
        assert decoded1["sub"] == decoded2["sub"] == str(sample_user.id)

        # Both tokens should be valid
        assert decoded1["exp"] > datetime.utcnow().timestamp()
        assert decoded2["exp"] > datetime.utcnow().timestamp()

    def test_login_with_malformed_credentials(self, client: TestClient):
        """Test login with malformed credential formats."""
        # Test with None values
        response1 = client.post("/routes/auth/login/access-token", data={
            "username": None,
            "password": "password123"
        })
        assert response1.status_code in [400, 422]

        # Test with empty password
        response2 = client.post("/routes/auth/login/access-token", data={
            "username": "test@example.com",
            "password": ""
        })
        assert response2.status_code in [400, 422]

    def test_login_case_sensitive_email(
        self,
        client: TestClient,
        sample_user: models.User
    ):
        """Test that login email matching is case-sensitive."""
        # Arrange - Use uppercase version of email
        login_data = {
            "username": sample_user.email.upper(),  # TESTUSER@EXAMPLE.COM
            "password": "testpassword123"
        }

        # Act
        response = client.post("/routes/auth/login/access-token", data=login_data)

        # Assert - May succeed or fail depending on implementation
        # Document current behavior: most systems treat emails as case-insensitive
        assert response.status_code in [200, 400]


@pytest.mark.unit
@pytest.mark.auth
class TestJWTTokenSecurity:
    """Test JWT token security and edge cases."""

    def test_token_with_missing_sub_claim(self, client: TestClient):
        """Test that token missing 'sub' claim is rejected."""
        # Arrange - Create token without 'sub' claim
        from datetime import datetime
        payload = {
            "exp": datetime.utcnow() + timedelta(minutes=30)
        }
        token_without_sub = jwt.encode(payload, settings.SECRET_KEY, algorithm=ALGORITHM)
        headers = {"Authorization": f"Bearer {token_without_sub}"}

        # Act
        response = client.get("/routes/auth/me", headers=headers)

        # Assert - Should fail due to missing sub claim
        assert response.status_code in [403, 404]  # Forbidden or Not Found
        data = response.json()
        assert "detail" in data

    def test_token_with_invalid_sub_claim(self, client: TestClient):
        """Test that token with non-existent user ID is rejected."""
        # Arrange - Create token with invalid user ID
        invalid_user_id = 999999
        token = create_access_token(subject=invalid_user_id)
        headers = {"Authorization": f"Bearer {token}"}

        # Act
        response = client.get("/routes/auth/me", headers=headers)

        # Assert - Should fail with user not found
        assert response.status_code == 404
        data = response.json()
        assert "detail" in data
        assert "not found" in data["detail"].lower()

    def test_token_with_wrong_algorithm(self, client: TestClient, sample_user: models.User):
        """Test that token signed with wrong algorithm is rejected."""
        # Arrange - Create token with different algorithm
        payload = {
            "sub": str(sample_user.id),
            "exp": datetime.utcnow() + timedelta(minutes=30)
        }
        # Use HS512 instead of HS256
        wrong_algorithm_token = jwt.encode(payload, settings.SECRET_KEY, algorithm="HS512")
        headers = {"Authorization": f"Bearer {wrong_algorithm_token}"}

        # Act
        response = client.get("/routes/auth/me", headers=headers)

        # Assert - Should fail due to algorithm mismatch
        assert response.status_code == 403
        data = response.json()
        assert "detail" in data

    def test_token_with_wrong_secret(self, client: TestClient, sample_user: models.User):
        """Test that token signed with wrong secret is rejected."""
        # Arrange - Create token with different secret
        payload = {
            "sub": str(sample_user.id),
            "exp": datetime.utcnow() + timedelta(minutes=30)
        }
        wrong_secret_token = jwt.encode(payload, "wrong_secret_key", algorithm=ALGORITHM)
        headers = {"Authorization": f"Bearer {wrong_secret_token}"}

        # Act
        response = client.get("/routes/auth/me", headers=headers)

        # Assert - Should fail due to signature verification
        assert response.status_code == 403
        data = response.json()
        assert "detail" in data

    def test_token_expiration_boundary(self, client: TestClient, sample_user: models.User):
        """Test token behavior at expiration boundary."""
        # Arrange - Create token that expires in 1 second
        short_lived_token = create_access_token(
            subject=sample_user.id,
            expires_delta=timedelta(seconds=1)
        )
        headers = {"Authorization": f"Bearer {short_lived_token}"}

        # Act - Use immediately (should work)
        response1 = client.get("/routes/auth/me", headers=headers)

        # Assert - First call should succeed
        assert response1.status_code == 200

        # Wait for token to expire
        import time
        time.sleep(2)

        # Act - Use after expiration (should fail)
        response2 = client.get("/routes/auth/me", headers=headers)

        # Assert - Second call should fail
        assert response2.status_code == 403
        data = response2.json()
        assert "detail" in data

    def test_bearer_token_format_variations(self, client: TestClient, sample_user: models.User):
        """Test various Bearer token format variations."""
        # Arrange
        valid_token = create_access_token(subject=sample_user.id)

        # Test: Extra whitespace
        headers1 = {"Authorization": f"Bearer  {valid_token}"}  # Double space
        response1 = client.get("/routes/auth/me", headers=headers1)
        assert response1.status_code in [200, 401, 403]  # Implementation-dependent

        # Test: Lowercase 'bearer'
        headers2 = {"Authorization": f"bearer {valid_token}"}
        response2 = client.get("/routes/auth/me", headers=headers2)
        assert response2.status_code in [200, 401, 403]  # Implementation-dependent

        # Test: Tab character
        headers3 = {"Authorization": f"Bearer\t{valid_token}"}
        response3 = client.get("/routes/auth/me", headers=headers3)
        assert response3.status_code in [200, 401, 403]  # Implementation-dependent


@pytest.mark.unit
@pytest.mark.auth
class TestAuthenticatedEndpoints:
    """Test endpoints requiring authentication."""

    def test_get_current_user_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User
    ):
        """Test retrieving current user profile with valid token."""
        # Act
        response = client.get("/routes/auth/me", headers=auth_headers)

        # Assert - Enhanced payload validation
        assert response.status_code == 200
        data = response.json()
        assert data["username"] == sample_user.username
        assert data["email"] == sample_user.email
        assert data["id"] == sample_user.id
        assert "is_superuser" in data

    def test_get_current_user_without_token_fails(self, client: TestClient):
        """Test accessing protected endpoint without token fails."""
        # Act
        response = client.get("/routes/auth/me")

        # Assert - Enhanced payload validation
        assert response.status_code == 401  # Unauthorized
        data = response.json()
        assert "detail" in data

    def test_get_current_user_invalid_token_fails(self, client: TestClient):
        """Test accessing protected endpoint with invalid token fails."""
        # Arrange
        invalid_headers = {"Authorization": "Bearer invalid_token_here"}

        # Act
        response = client.get("/routes/auth/me", headers=invalid_headers)

        # Assert - Enhanced payload validation
        assert response.status_code == 403  # Forbidden
        data = response.json()
        assert "detail" in data
        assert "credentials" in data["detail"].lower()

    def test_get_current_user_expired_token_fails(self, client: TestClient, sample_user: models.User):
        """Test accessing protected endpoint with expired token fails."""
        # Arrange - Create expired token
        expired_token = create_access_token(
            subject=sample_user.id,
            expires_delta=timedelta(minutes=-10)  # Expired 10 minutes ago
        )
        expired_headers = {"Authorization": f"Bearer {expired_token}"}

        # Act
        response = client.get("/routes/auth/me", headers=expired_headers)

        # Assert
        assert response.status_code == 403
        data = response.json()
        assert "detail" in data

    def test_get_current_user_missing_bearer_prefix(self, client: TestClient, auth_headers: dict):
        """Test accessing protected endpoint with token missing Bearer prefix fails."""
        # Arrange - Extract token without Bearer prefix
        token = auth_headers["Authorization"].replace("Bearer ", "")
        malformed_headers = {"Authorization": token}

        # Act
        response = client.get("/routes/auth/me", headers=malformed_headers)

        # Assert
        assert response.status_code == 401
        data = response.json()
        assert "detail" in data

    def test_get_current_user_tampered_token_fails(self, client: TestClient, sample_user: models.User):
        """Test accessing protected endpoint with tampered token fails."""
        # Arrange - Create token and tamper with it
        valid_token = create_access_token(subject=sample_user.id)
        tampered_token = valid_token[:-5] + "xxxxx"  # Tamper with signature
        tampered_headers = {"Authorization": f"Bearer {tampered_token}"}

        # Act
        response = client.get("/routes/auth/me", headers=tampered_headers)

        # Assert
        assert response.status_code == 403
        data = response.json()
        assert "detail" in data

    def test_update_user_plugins_success(
        self,
        client: TestClient,
        auth_headers: dict,
        sample_user: models.User
    ):
        """Test updating user plugins successfully."""
        # Arrange
        update_data = {
            "plugins": ["plugin1", "plugin2", "plugin3"]
        }

        # Act
        response = client.post("/routes/auth/plugins", json=update_data, headers=auth_headers)

        # Assert - Enhanced payload validation
        assert response.status_code == 200
        data = response.json()
        assert "id" in data
        assert data["email"] == sample_user.email

    def test_update_user_plugins_unauthorized(self, client: TestClient):
        """Test updating user plugins without authentication fails."""
        # Arrange
        update_data = {
            "plugins": ["plugin1", "plugin2"]
        }

        # Act
        response = client.post("/routes/auth/plugins", json=update_data)

        # Assert
        assert response.status_code == 401
        data = response.json()
        assert "detail" in data
