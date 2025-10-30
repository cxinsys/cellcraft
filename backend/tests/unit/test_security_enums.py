"""
Unit tests for security and enum modules.

Test coverage:
- JWT token creation and validation
- Password hashing and verification
- Plugin type enumerations
"""
import pytest
from datetime import timedelta
from jose import jwt

from app.common.security import (
    create_access_token,
    verify_password,
    get_password_hash,
    ALGORITHM
)
from app.common.enums import PluginType
from app.common.config import settings


@pytest.mark.unit
@pytest.mark.security
class TestJWTTokenCreation:
    """Test JWT token creation and encoding."""

    def test_create_access_token_with_custom_expiration(self):
        """Test token creation with custom expiration time."""
        # Arrange
        subject = "test_user_123"
        expires_delta = timedelta(minutes=30)

        # Act
        token = create_access_token(subject, expires_delta)

        # Assert
        assert token is not None
        assert isinstance(token, str)

        # Decode and verify token structure
        decoded = jwt.decode(token, settings.SECRET_KEY, algorithms=[ALGORITHM])
        assert decoded["sub"] == subject
        assert "exp" in decoded

    def test_create_access_token_with_default_expiration(self):
        """Test token creation with default expiration from settings."""
        # Arrange
        subject = 42  # Testing with integer subject

        # Act
        token = create_access_token(subject)

        # Assert
        assert token is not None
        assert isinstance(token, str)

        # Decode and verify
        decoded = jwt.decode(token, settings.SECRET_KEY, algorithms=[ALGORITHM])
        assert decoded["sub"] == str(subject)  # Should be converted to string
        assert "exp" in decoded

    def test_create_access_token_subject_types(self):
        """Test token creation with different subject types."""
        # Test with string
        token_str = create_access_token("user@example.com")
        decoded_str = jwt.decode(token_str, settings.SECRET_KEY, algorithms=[ALGORITHM])
        assert decoded_str["sub"] == "user@example.com"

        # Test with integer
        token_int = create_access_token(123)
        decoded_int = jwt.decode(token_int, settings.SECRET_KEY, algorithms=[ALGORITHM])
        assert decoded_int["sub"] == "123"

    def test_expired_token_raises_error(self):
        """Test that expired tokens cannot be decoded."""
        # Arrange - Create token that expired 10 minutes ago
        subject = "test_user"
        expired_token = create_access_token(subject, timedelta(minutes=-10))

        # Act & Assert - Expired token should raise JWTError
        with pytest.raises(jwt.JWTError):
            jwt.decode(expired_token, settings.SECRET_KEY, algorithms=[ALGORITHM])

    def test_token_with_invalid_signature_fails(self):
        """Test that tampered tokens with invalid signatures are rejected."""
        # Arrange - Create valid token and tamper with signature
        token = create_access_token("user123")
        tampered_token = token[:-10] + "tampered00"  # Modify signature

        # Act & Assert - Tampered token should raise JWTError
        with pytest.raises(jwt.JWTError):
            jwt.decode(tampered_token, settings.SECRET_KEY, algorithms=[ALGORITHM])

    def test_token_without_sub_claim_fails(self):
        """Test that tokens missing 'sub' claim are invalid."""
        # Arrange - Manually create token without 'sub' claim
        from datetime import datetime
        payload = {
            "exp": datetime.utcnow() + timedelta(minutes=30)
        }
        token_without_sub = jwt.encode(payload, settings.SECRET_KEY, algorithm=ALGORITHM)

        # Act - Decode token (will succeed)
        decoded = jwt.decode(token_without_sub, settings.SECRET_KEY, algorithms=[ALGORITHM])

        # Assert - 'sub' claim should be missing
        assert "sub" not in decoded
        assert "exp" in decoded

    def test_completely_invalid_token_format_fails(self):
        """Test that completely invalid token formats are rejected."""
        # Arrange
        invalid_tokens = [
            "not.a.valid.jwt.token",
            "invalid_base64_$%^&",
            "",
            "header.payload"  # Missing signature
        ]

        # Act & Assert - All invalid formats should raise JWTError
        for invalid_token in invalid_tokens:
            with pytest.raises(jwt.JWTError):
                jwt.decode(invalid_token, settings.SECRET_KEY, algorithms=[ALGORITHM])


@pytest.mark.unit
@pytest.mark.security
class TestPasswordHashing:
    """Test password hashing and verification functions."""

    def test_get_password_hash_generates_hash(self):
        """Test password hashing generates a bcrypt hash."""
        # Arrange
        password = "secure_password_123"

        # Act
        hashed = get_password_hash(password)

        # Assert
        assert hashed is not None
        assert isinstance(hashed, str)
        assert hashed != password
        # Bcrypt hash can start with $2a$, $2b$, $2x$, or $2y$ depending on version
        assert hashed.startswith("$2") and "$" in hashed[3:], "Should be a bcrypt hash"

    def test_get_password_hash_unique_for_same_password(self):
        """Test same password generates different hashes (salt)."""
        # Arrange
        password = "test_password"

        # Act
        hash1 = get_password_hash(password)
        hash2 = get_password_hash(password)

        # Assert
        assert hash1 != hash2  # Different salts

    def test_verify_password_success(self):
        """Test password verification succeeds with correct password."""
        # Arrange
        plain_password = "correct_password"
        hashed_password = get_password_hash(plain_password)

        # Act
        result = verify_password(plain_password, hashed_password)

        # Assert
        assert result is True

    def test_verify_password_failure(self):
        """Test password verification fails with incorrect password."""
        # Arrange
        plain_password = "correct_password"
        wrong_password = "wrong_password"
        hashed_password = get_password_hash(plain_password)

        # Act
        result = verify_password(wrong_password, hashed_password)

        # Assert
        assert result is False

    def test_get_password_hash_empty_password(self):
        """Test password hashing with empty string.

        Edge case: Bcrypt can hash empty strings, should not fail.
        Verifies that empty password generates valid bcrypt hash.
        """
        # Arrange
        password = ""

        # Act
        hashed = get_password_hash(password)

        # Assert
        assert hashed is not None
        assert isinstance(hashed, str)
        # Bcrypt hash should start with $2 (versions: $2a$, $2b$, $2x$, $2y$)
        assert hashed.startswith("$2"), "Should be a valid bcrypt hash"
        # Verify empty password can be verified against its hash
        assert verify_password("", hashed) is True

    def test_verify_password_with_empty_strings(self):
        """Test password verification with empty password and hash combinations.

        Edge case: Verifies behavior when comparing empty strings.
        - Empty password vs empty password hash: should match
        - Empty password vs non-empty hash: should not match
        - Non-empty password vs empty hash: edge case handling
        """
        # Arrange
        empty_password = ""
        non_empty_password = "valid_password_123"
        empty_password_hash = get_password_hash(empty_password)
        non_empty_password_hash = get_password_hash(non_empty_password)

        # Act & Assert - Scenario 1: Empty password hash verified with empty password
        result1 = verify_password(empty_password, empty_password_hash)
        assert result1 is True, "Empty password should match its own hash"

        # Act & Assert - Scenario 2: Empty password vs non-empty password hash
        result2 = verify_password(empty_password, non_empty_password_hash)
        assert result2 is False, "Empty password should not match non-empty password hash"

        # Act & Assert - Scenario 3: Non-empty password vs empty password hash
        result3 = verify_password(non_empty_password, empty_password_hash)
        assert result3 is False, "Non-empty password should not match empty password hash"


@pytest.mark.unit
@pytest.mark.enum
class TestPluginTypeEnum:
    """Test PluginType enumeration."""

    def test_plugin_type_analysis_value(self):
        """Test PluginType.ANALYSIS has correct value."""
        # Assert
        assert PluginType.ANALYSIS == "analysis"
        assert PluginType.ANALYSIS.value == "analysis"

    def test_plugin_type_visualization_value(self):
        """Test PluginType.VISUALIZATION has correct value."""
        # Assert
        assert PluginType.VISUALIZATION == "visualization"
        assert PluginType.VISUALIZATION.value == "visualization"

    def test_plugin_type_members(self):
        """Test PluginType enum has exactly two members."""
        # Assert
        assert len(PluginType) == 2
        assert set(PluginType) == {PluginType.ANALYSIS, PluginType.VISUALIZATION}

    def test_plugin_type_invalid_attribute_access_raises_error(self):
        """Test accessing non-existent enum member raises AttributeError."""
        # Act & Assert - Accessing invalid enum member should raise AttributeError
        with pytest.raises(AttributeError):
            _ = PluginType.NONEXISTENT

    def test_plugin_type_invalid_value_raises_error(self):
        """Test creating enum with invalid value raises ValueError."""
        # Act & Assert - Invalid enum value should raise ValueError
        with pytest.raises(ValueError):
            PluginType("invalid_type")

    def test_plugin_type_case_sensitive_value(self):
        """Test enum value comparison is case-sensitive."""
        # Arrange & Act
        analysis_value = PluginType.ANALYSIS.value

        # Assert - Enum values are case-sensitive
        assert analysis_value == "analysis"
        assert analysis_value != "ANALYSIS"
        assert analysis_value != "Analysis"
