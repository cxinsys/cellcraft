"""
Unit tests for Pydantic schema validation.

Test coverage:
- UserCreate schema validation (required fields, types, BUG-001, BUG-002, BUG-003)
- UserUpdate schema validation (optional fields, partial updates)
- UserProfile schema validation (response serialization)
- UserBase schema validation (shared properties)

Related bugs:
- BUG-001: None values bypass Pydantic validation
- BUG-002: Empty string validation gap for username
- BUG-003: Username field made required in UserCreate schema
"""
import pytest
from pydantic import ValidationError

from app.user.schemas import (
    UserBase,
    UserCreate,
    UserUpdate,
    UserProfile,
    User
)


@pytest.mark.unit
@pytest.mark.api
class TestUserCreateSchema:
    """Test UserCreate schema validation for user registration.

    UserCreate is used for POST /routes/auth/register endpoint.
    Critical for data integrity and API security.

    Required fields (BUG-003 fix): username, email, password
    """

    def test_user_create_valid_data(self):
        """Test UserCreate accepts valid registration data."""
        # Arrange
        valid_data = {
            "username": "testuser",
            "email": "test@example.com",
            "password": "securepassword123"
        }

        # Act
        user = UserCreate(**valid_data)

        # Assert
        assert user.username == "testuser"
        assert user.email == "test@example.com"
        assert user.password == "securepassword123"
        assert user.is_active is True  # Default value
        assert user.is_superuser is False  # Default value

    def test_user_create_missing_username_raises_validation_error(self):
        """Test UserCreate rejects missing username field.

        Previously username was optional, but BUG-003 fix made it required.
        Expected: ValidationError (422 Unprocessable Entity)
        """
        # Arrange
        invalid_data = {
            "email": "test@example.com",
            "password": "securepassword123"
            # Missing required username field
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("username",) for error in errors)
        assert any("field required" in error["msg"].lower() or
                   "missing" in error["msg"].lower()
                   for error in errors)

    def test_user_create_none_email_raises_validation_error(self):
        """Test UserCreate rejects None for required email field.

        Previously BUG-001: Fixed by Pydantic v2 - None values are now caught.
        Expected: ValidationError (422 Unprocessable Entity)
        """
        # Arrange
        invalid_data = {
            "username": "testuser",
            "email": None,  # Required field
            "password": "securepassword123"
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("email",) for error in errors)
        assert any("none is not an allowed value" in error["msg"].lower() or
                   "field cannot be none" in error["msg"].lower()
                   for error in errors)

    def test_user_create_none_password_raises_validation_error(self):
        """Test UserCreate rejects None for required password field.

        Previously BUG-001: Fixed by Pydantic v2 - None values are now caught.
        """
        # Arrange
        invalid_data = {
            "username": "testuser",
            "email": "test@example.com",
            "password": None  # Required field
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("password",) for error in errors)

    def test_user_create_empty_username_raises_validation_error(self):
        """Test UserCreate rejects empty string for username field.

        Previously BUG-002: Fixed by adding Field(min_length=1) and custom validator.
        Expected: ValidationError (min_length constraint)
        """
        # Arrange
        invalid_data = {
            "username": "",  # Empty string should be rejected
            "email": "test@example.com",
            "password": "securepassword123"
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("username",) for error in errors)
        assert any("cannot be empty" in error["msg"].lower() or
                   "at least" in error["msg"].lower() or
                   "min_length" in error["msg"].lower()
                   for error in errors)

    def test_user_create_empty_password_raises_validation_error(self):
        """Test UserCreate rejects empty string for password field.

        Previously BUG-002: Fixed by adding Field(min_length=1) and custom validator.
        """
        # Arrange
        invalid_data = {
            "username": "testuser",
            "email": "test@example.com",
            "password": ""  # Empty password
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("password",) for error in errors)

    def test_user_create_invalid_email_format(self):
        """Test UserCreate rejects invalid email format."""
        # Arrange
        invalid_data = {
            "username": "testuser",
            "email": "not-an-email",  # Invalid format
            "password": "securepassword123"
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("email",) for error in errors)
        assert any("email" in error["msg"].lower() for error in errors)

    def test_user_create_email_with_plus_sign(self):
        """Test UserCreate accepts email with plus sign (valid RFC 5322)."""
        # Arrange
        valid_data = {
            "username": "testuser",
            "email": "test+tag@example.com",  # Plus sign is valid
            "password": "securepassword123"
        }

        # Act
        user = UserCreate(**valid_data)

        # Assert
        assert user.email == "test+tag@example.com"

    def test_user_create_missing_required_fields(self):
        """Test UserCreate raises ValidationError when required fields are missing."""
        # Arrange - missing email and password
        invalid_data = {
            "username": "testuser"
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserCreate(**invalid_data)

        # Verify multiple fields are reported as missing
        errors = exc_info.value.errors()
        missing_fields = [error["loc"][0] for error in errors]
        assert "email" in missing_fields
        assert "password" in missing_fields

    def test_user_create_extra_fields_ignored(self):
        """Test UserCreate ignores extra fields not in schema."""
        # Arrange
        data_with_extra = {
            "username": "testuser",
            "email": "test@example.com",
            "password": "securepassword123",
            "extra_field": "should_be_ignored"
        }

        # Act
        user = UserCreate(**data_with_extra)

        # Assert - extra field should not exist
        assert not hasattr(user, "extra_field")
        assert user.username == "testuser"


@pytest.mark.unit
@pytest.mark.api
class TestUserUpdateSchema:
    """Test UserUpdate schema validation for user profile updates.

    UserUpdate is used for PATCH /routes/users/{user_id} endpoint.
    All fields are optional for partial updates.
    """

    def test_user_update_all_fields(self):
        """Test UserUpdate accepts all optional fields."""
        # Arrange
        update_data = {
            "username": "newusername",
            "email": "newemail@example.com",
            "password": "newpassword123",
            "is_active": False,
            "is_superuser": True,
            "plugins": [1, 2, 3]
        }

        # Act
        user_update = UserUpdate(**update_data)

        # Assert
        assert user_update.username == "newusername"
        assert user_update.email == "newemail@example.com"
        assert user_update.password == "newpassword123"
        assert user_update.is_active is False
        assert user_update.is_superuser is True
        assert user_update.plugins == [1, 2, 3]

    def test_user_update_partial_fields(self):
        """Test UserUpdate accepts partial updates (only some fields)."""
        # Arrange - only update email
        update_data = {
            "email": "newemail@example.com"
        }

        # Act
        user_update = UserUpdate(**update_data)

        # Assert
        assert user_update.email == "newemail@example.com"
        assert user_update.username is None
        assert user_update.password is None
        assert user_update.is_active is True  # Default from UserBase

    def test_user_update_empty_dict(self):
        """Test UserUpdate accepts empty dict (no fields to update)."""
        # Arrange
        update_data = {}

        # Act
        user_update = UserUpdate(**update_data)

        # Assert - all fields should be None or defaults
        assert user_update.username is None
        assert user_update.email is None
        assert user_update.password is None
        assert user_update.plugins is None

    def test_user_update_invalid_email_format(self):
        """Test UserUpdate rejects invalid email format."""
        # Arrange
        invalid_data = {
            "email": "not-an-email"
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserUpdate(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("email",) for error in errors)

    def test_user_update_set_field_to_none(self):
        """Test UserUpdate allows setting optional fields to None."""
        # Arrange
        update_data = {
            "username": None,  # Explicitly set to None
            "email": "test@example.com"
        }

        # Act
        user_update = UserUpdate(**update_data)

        # Assert
        assert user_update.username is None
        assert user_update.email == "test@example.com"

    def test_user_update_plugins_list_validation(self):
        """Test UserUpdate validates plugins field as list."""
        # Arrange
        update_data = {
            "plugins": [1, 2, 3, 4]
        }

        # Act
        user_update = UserUpdate(**update_data)

        # Assert
        assert isinstance(user_update.plugins, list)
        assert len(user_update.plugins) == 4


@pytest.mark.unit
@pytest.mark.api
class TestUserProfileSchema:
    """Test UserProfile schema validation for API responses.

    UserProfile is used for GET /routes/users/me endpoint.
    Should include all public user information.
    """

    def test_user_profile_valid_data(self):
        """Test UserProfile serializes valid user data."""
        # Arrange
        profile_data = {
            "id": 1,
            "username": "testuser",
            "email": "test@example.com",
            "is_active": True,
            "is_superuser": False
        }

        # Act
        profile = UserProfile(**profile_data)

        # Assert
        assert profile.id == 1
        assert profile.username == "testuser"
        assert profile.email == "test@example.com"
        assert profile.is_active is True
        assert profile.is_superuser is False

    def test_user_profile_missing_required_fields(self):
        """Test UserProfile requires id, username, email, is_superuser."""
        # Arrange - missing email
        invalid_data = {
            "id": 1,
            "username": "testuser",
            "is_superuser": False
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            UserProfile(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("email",) for error in errors)

    def test_user_profile_orm_mode_enabled(self):
        """Test UserProfile Config has orm_mode enabled for ORM compatibility."""
        # Assert
        assert hasattr(UserProfile, 'Config')
        assert hasattr(UserProfile.Config, 'orm_mode')
        assert UserProfile.Config.orm_mode is True

    def test_user_profile_excludes_password(self):
        """Test UserProfile does not include password or hashed_password fields."""
        # Arrange
        profile_data = {
            "id": 1,
            "username": "testuser",
            "email": "test@example.com",
            "is_active": True,
            "is_superuser": False,
            "password": "should_not_appear",  # Should be ignored
            "hashed_password": "should_not_appear"  # Should be ignored
        }

        # Act
        profile = UserProfile(**profile_data)

        # Assert - password fields should not exist
        assert not hasattr(profile, "password")
        assert not hasattr(profile, "hashed_password")


@pytest.mark.unit
@pytest.mark.api
class TestUserSchema:
    """Test User schema validation (internal representation with password).

    User schema includes hashed_password for internal use.
    Should not be exposed in API responses.
    """

    def test_user_schema_valid_data(self):
        """Test User schema accepts valid data with hashed_password."""
        # Arrange
        user_data = {
            "id": 1,
            "username": "testuser",
            "email": "test@example.com",
            "is_active": True,
            "is_superuser": False,
            "hashed_password": "$2b$12$hashed_password_string"
        }

        # Act
        user = User(**user_data)

        # Assert
        assert user.id == 1
        assert user.username == "testuser"
        assert user.email == "test@example.com"
        assert user.is_active is True
        assert user.is_superuser is False
        assert user.hashed_password == "$2b$12$hashed_password_string"

    def test_user_schema_requires_hashed_password(self):
        """Test User schema requires hashed_password field."""
        # Arrange - missing hashed_password
        invalid_data = {
            "id": 1,
            "username": "testuser",
            "email": "test@example.com",
            "is_active": True,
            "is_superuser": False
        }

        # Act & Assert
        with pytest.raises(ValidationError) as exc_info:
            User(**invalid_data)

        # Verify error details
        errors = exc_info.value.errors()
        assert any(error["loc"] == ("hashed_password",) for error in errors)

    def test_user_schema_orm_mode_enabled(self):
        """Test User schema Config has orm_mode enabled for ORM compatibility."""
        # Assert
        assert hasattr(User, 'Config')
        assert hasattr(User.Config, 'orm_mode')
        assert User.Config.orm_mode is True


@pytest.mark.unit
@pytest.mark.api
class TestUserBaseSchema:
    """Test UserBase shared properties schema.

    UserBase is inherited by other user schemas.
    Defines common optional fields with defaults.
    """

    def test_user_base_default_values(self):
        """Test UserBase provides correct default values."""
        # Arrange & Act
        base = UserBase()

        # Assert
        assert base.email is None
        assert base.username is None
        assert base.is_active is True  # Default True
        assert base.is_superuser is False  # Default False

    def test_user_base_accepts_all_fields(self):
        """Test UserBase accepts all defined fields."""
        # Arrange
        base_data = {
            "email": "test@example.com",
            "username": "testuser",
            "is_active": False,
            "is_superuser": True
        }

        # Act
        base = UserBase(**base_data)

        # Assert
        assert base.email == "test@example.com"
        assert base.username == "testuser"
        assert base.is_active is False
        assert base.is_superuser is True

    def test_user_base_accepts_partial_fields(self):
        """Test UserBase accepts partial field updates."""
        # Arrange
        base_data = {
            "email": "test@example.com"
        }

        # Act
        base = UserBase(**base_data)

        # Assert
        assert base.email == "test@example.com"
        assert base.username is None
        assert base.is_active is True  # Default
        assert base.is_superuser is False  # Default
