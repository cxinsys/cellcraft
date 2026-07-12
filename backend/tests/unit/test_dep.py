"""
Unit tests for dependency injection functions.

Test coverage:
- Database session management (get_db)
- User authentication and token validation (get_current_user)
- Active user validation (get_current_active_user)
"""
import pytest
from unittest.mock import MagicMock, patch
from datetime import datetime, timedelta
from fastapi import HTTPException
from jose import jwt

from app.db.session import get_db
from app.routes.dep import get_current_user, get_current_active_user
from app.database import models
from app.core.config import settings
from app.core.security import ALGORITHM, create_access_token


@pytest.mark.unit
@pytest.mark.auth
class TestGetDb:
    """Test database session dependency."""

    @patch('app.db.session.SessionLocal')
    def test_get_db_yields_session(self, mock_session_local):
        """Test get_db yields database session and closes it."""
        # Arrange
        mock_db = MagicMock()
        mock_session_local.return_value = mock_db

        # Act
        db_generator = get_db()
        db = next(db_generator)

        # Assert
        assert db == mock_db
        mock_session_local.assert_called_once()

        # Test cleanup - close is called when generator exits
        try:
            next(db_generator)
        except StopIteration:
            pass

        mock_db.close.assert_called_once()

    @patch('app.db.session.SessionLocal')
    def test_get_db_closes_on_exception(self, mock_session_local):
        """Test get_db closes session even if exception occurs."""
        # Arrange
        mock_db = MagicMock()
        mock_session_local.return_value = mock_db

        # Act
        db_generator = get_db()
        db = next(db_generator)

        # Simulate exception during usage
        try:
            db_generator.throw(Exception("Test exception"))
        except Exception:
            pass

        # Assert - close should still be called
        mock_db.close.assert_called_once()


@pytest.mark.unit
@pytest.mark.auth
class TestGetCurrentUser:
    """Test current user authentication dependency."""

    def test_get_current_user_success(self, monkeypatch):
        """Test successful user retrieval with valid token using isolated SECRET_KEY."""
        # Arrange - Use monkeypatch for SECRET_KEY isolation
        test_secret_key = "test_secret_key_for_isolation_12345678"
        monkeypatch.setattr('app.routes.dep.settings.SECRET_KEY', test_secret_key)

        mock_db = MagicMock()
        user_id = 123
        token = jwt.encode({"sub": str(user_id)}, test_secret_key, algorithm=ALGORITHM)

        mock_user = models.User(
            id=user_id,
            email="test@example.com",
            username="testuser",
            hashed_password="hashed",
            is_active=True,
            is_superuser=False
        )

        with patch('app.routes.dep.crud_user.get_user', return_value=mock_user):
            # Act
            result = get_current_user(db=mock_db, token=token)

            # Assert
            assert result == mock_user
            assert result.id == user_id

    def test_get_current_user_invalid_token(self):
        """Test get_current_user raises 403 for invalid token."""
        # Arrange
        mock_db = MagicMock()
        invalid_token = "invalid.jwt.token"

        # Act & Assert
        with pytest.raises(HTTPException) as exc_info:
            get_current_user(db=mock_db, token=invalid_token)

        assert exc_info.value.status_code == 403
        assert exc_info.value.detail == "Could not validate credentials"

    def test_get_current_user_expired_token(self):
        """Test get_current_user raises 403 for actually expired token."""
        # Arrange
        mock_db = MagicMock()
        user_id = 123

        # Create a token that expired 10 minutes ago using real expiration
        expired_token = create_access_token(subject=user_id, expires_delta=timedelta(minutes=-10))

        # Act & Assert - Expired token should raise 403
        with pytest.raises(HTTPException) as exc_info:
            get_current_user(db=mock_db, token=expired_token)

        assert exc_info.value.status_code == 403
        assert exc_info.value.detail == "Could not validate credentials"

    def test_get_current_user_token_without_sub_claim(self):
        """Test get_current_user behavior when token missing 'sub' claim.

        Note: This test documents current behavior where missing 'sub' results in
        user_id=None, which then causes 404 when user lookup fails.
        A future improvement would be to validate 'sub' claim existence early.
        """
        # Arrange
        mock_db = MagicMock()

        # Manually create token without 'sub' claim
        payload = {"exp": datetime.utcnow() + timedelta(minutes=30)}
        token_without_sub = jwt.encode(payload, settings.SECRET_KEY, algorithm=ALGORITHM)

        # Mock the user lookup to return None (user not found)
        with patch('app.routes.dep.crud_user.get_user', return_value=None):
            # Act & Assert - Missing sub leads to None user_id, then 404
            with pytest.raises(HTTPException) as exc_info:
                get_current_user(db=mock_db, token=token_without_sub)

            assert exc_info.value.status_code == 404
            assert exc_info.value.detail == "User not found"

    def test_get_current_user_not_found(self):
        """Test get_current_user raises 404 when user doesn't exist."""
        # Arrange
        mock_db = MagicMock()
        user_id = 999
        token = jwt.encode({"sub": str(user_id)}, settings.SECRET_KEY, algorithm=ALGORITHM)

        with patch('app.routes.dep.crud_user.get_user', return_value=None):
            # Act & Assert
            with pytest.raises(HTTPException) as exc_info:
                get_current_user(db=mock_db, token=token)

            assert exc_info.value.status_code == 404
            assert exc_info.value.detail == "User not found"

    def test_get_current_user_malformed_token_structure(self):
        """Test get_current_user with malformed JWT token structures.

        JWT tokens should have 3 parts: header.payload.signature
        This test validates proper error handling for:
        - Token with only 2 parts (missing signature)
        - Token with 4+ parts (extra segments)
        - Token with special characters in wrong places

        Expected: All should raise HTTPException with 403 status
        """
        # Arrange
        mock_db = MagicMock()
        malformed_tokens = [
            "header.payload",  # Only 2 parts (missing signature)
            "header.payload.signature.extra",  # 4 parts (extra segments)
            "not-a-valid-jwt"  # Invalid format
        ]

        for token in malformed_tokens:
            # Act & Assert
            with pytest.raises(HTTPException) as exc_info:
                get_current_user(db=mock_db, token=token)

            assert exc_info.value.status_code == 403
            assert exc_info.value.detail == "Could not validate credentials"

@pytest.mark.unit
@pytest.mark.auth
class TestGetCurrentActiveUser:
    """Test active user validation dependency."""

    def test_get_current_active_user_success(self):
        """Test successful retrieval of active user."""
        # Arrange
        mock_user = models.User(
            id=123,
            email="active@example.com",
            username="activeuser",
            hashed_password="hashed",
            is_active=True,
            is_superuser=False
        )

        with patch('app.routes.dep.crud_user.is_active', return_value=True):
            # Act
            result = get_current_active_user(current_user=mock_user)

            # Assert
            assert result == mock_user
            assert result.is_active is True

    def test_get_current_active_user_inactive(self):
        """Test get_current_active_user raises 400 for inactive user."""
        # Arrange
        mock_user = models.User(
            id=123,
            email="inactive@example.com",
            username="inactiveuser",
            hashed_password="hashed",
            is_active=False,
            is_superuser=False
        )

        with patch('app.routes.dep.crud_user.is_active', return_value=False):
            # Act & Assert
            with pytest.raises(HTTPException) as exc_info:
                get_current_active_user(current_user=mock_user)

            assert exc_info.value.status_code == 400
            assert exc_info.value.detail == "Inactive user"

    def test_get_current_active_user_error_message_consistency(self):
        """Test error message consistency for inactive user validation.

        Ensures that the error detail message is exactly "Inactive user"
        for proper frontend error handling and user messaging.

        This validates the exact string used in dep.py:53
        """
        # Arrange
        mock_user = models.User(
            id=456,
            email="consistency@example.com",
            username="consistencyuser",
            hashed_password="hashed",
            is_active=False,
            is_superuser=False
        )

        with patch('app.routes.dep.crud_user.is_active', return_value=False):
            # Act & Assert
            with pytest.raises(HTTPException) as exc_info:
                get_current_active_user(current_user=mock_user)

            # Assert - Verify exact error message consistency
            assert exc_info.value.status_code == 400
            assert exc_info.value.detail == "Inactive user"
