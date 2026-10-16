"""
Unit tests for CRUD user operations (app/user/crud.py).

Test coverage:
- Create operations: create_user (with plugin auto-association), create_superuser
- Read operations: get_user, get_user_by_email, get_users_count
- Update operations: update_user (password updates)
- Authentication: authenticate, is_active, is_superuser

Total: 24 tests across 5 test classes

Critical Business Logic Tested:
- Plugin auto-association when creating users (all plugins automatically linked)
- Password hashing verification
- Authentication credential verification
- User status checks (active/superuser)
"""
import pytest
from sqlalchemy.orm import Session

from app.user import crud as crud_user
from app import models
from app.user.schemas import UserCreate, UserUpdate
from app.core.security import verify_password


@pytest.mark.unit
@pytest.mark.crud
class TestCreateUser:
    """Test create_user() function with plugin auto-association logic.

    Critical: create_user() automatically associates ALL plugins to new users.
    This must be verified to prevent user-plugin desynchronization.
    """

    def test_create_user_success(self, db_session: Session):
        """Test successful user creation with all required fields."""
        # Arrange
        user_data = UserCreate(
            username="newuser",
            email="newuser@example.com",
            password="securepassword123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)

        # Assert
        assert created_user.id is not None
        assert created_user.username == "newuser"
        assert created_user.email == "newuser@example.com"
        assert created_user.hashed_password is not None
        assert created_user.is_active is True
        assert created_user.is_superuser is False

    def test_create_user_with_plugins_auto_association(self, db_session: Session, sample_plugin):
        """Test that create_user() auto-associates ALL plugins to new user.

        Business Rule: When a user is created, they should automatically
        have access to all plugins in the system.
        """
        # Arrange - sample_plugin already exists in DB
        user_data = UserCreate(
            username="testuser",
            email="testuser@example.com",
            password="password123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)
        db_session.refresh(created_user)

        # Assert - User should have the plugin associated
        assert len(created_user.plugins) == 1
        assert created_user.plugins[0].id == sample_plugin.id
        assert created_user.plugins[0].name == sample_plugin.name

    def test_create_user_with_zero_plugins(self, db_session: Session):
        """Test user creation when no plugins exist in database (edge case)."""
        # Arrange - No plugins in DB
        user_data = UserCreate(
            username="testuser",
            email="testuser@example.com",
            password="password123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)
        db_session.refresh(created_user)

        # Assert - User created successfully with empty plugins list
        assert created_user.id is not None
        assert len(created_user.plugins) == 0

    def test_create_user_with_multiple_plugins(self, db_session: Session):
        """Test plugin auto-association with 5+ plugins."""
        # Arrange - Create 5 plugins
        from app.core.enums import PluginType

        for i in range(5):
            plugin = models.Plugin(
                name=f"Plugin{i+1}",
                description=f"Test plugin {i+1}",
                author="test",
                plugin_path=f"./plugin/local/Plugin{i+1}/",
                plugin_type=PluginType.ANALYSIS,
                dependencies={},
                drawflow={},
                rules={},
                use_gpu=False,
                source="local",
                version="1.0.0"
            )
            db_session.add(plugin)
        db_session.commit()

        user_data = UserCreate(
            username="testuser",
            email="testuser@example.com",
            password="password123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)
        db_session.refresh(created_user)

        # Assert - All 5 plugins associated
        assert len(created_user.plugins) == 5
        plugin_names = [p.name for p in created_user.plugins]
        assert set(plugin_names) == {"Plugin1", "Plugin2", "Plugin3", "Plugin4", "Plugin5"}

    def test_create_user_hashes_password(self, db_session: Session):
        """Test that password is hashed and not stored in plaintext."""
        # Arrange
        plain_password = "my_secret_password_123"
        user_data = UserCreate(
            username="secureuser",
            email="secure@example.com",
            password=plain_password
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)

        # Assert - Password should be hashed
        assert created_user.hashed_password != plain_password
        assert created_user.hashed_password.startswith("$2b$")  # bcrypt hash format
        assert len(created_user.hashed_password) == 60  # bcrypt hash length
        # Verify password can be verified with hash
        assert verify_password(plain_password, created_user.hashed_password) is True

    def test_create_user_with_explicit_username(self, db_session: Session):
        """Test user creation requires username (database and schema constraint).

        Schema-Model Consistency:
        - Schema (schemas/user.py:15): username: str = Field(..., min_length=1) [FIXED]
        - Model (models.py:19): username = Column(String, nullable=False)

        Both schema and model now require username.
        """
        # Arrange
        user_data = UserCreate(
            username="requireduser",  # Required by both schema and database
            email="required@example.com",
            password="password123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)

        # Assert
        assert created_user.id is not None
        assert created_user.username == "requireduser"
        assert created_user.email == "required@example.com"

    def test_create_user_default_values(self, db_session: Session):
        """Test that new users have correct default values for is_active and is_superuser."""
        # Arrange
        user_data = UserCreate(
            username="defaultuser",
            email="default@example.com",
            password="password123"
        )

        # Act
        created_user = crud_user.create_user(db_session, user=user_data)

        # Assert - Verify defaults from User model
        assert created_user.is_active is True  # Should be active by default
        assert created_user.is_superuser is False  # Should not be superuser


@pytest.mark.unit
@pytest.mark.crud
class TestCreateSuperuser:
    """Test create_superuser() function for admin user creation."""

    def test_create_superuser_success(self, db_session: Session):
        """Test successful superuser creation with elevated privileges."""
        # Arrange
        admin_data = UserCreate(
            username="admin",
            email="admin@example.com",
            password="adminpassword123"
        )

        # Act
        created_admin = crud_user.create_superuser(db_session, user=admin_data)

        # Assert
        assert created_admin.id is not None
        assert created_admin.username == "admin"
        assert created_admin.email == "admin@example.com"
        assert created_admin.is_superuser is True  # Critical: must be superuser
        assert created_admin.is_active is True

    def test_create_superuser_is_active_default(self, db_session: Session):
        """Test that superusers are always active by default."""
        # Arrange
        admin_data = UserCreate(
            username="superadmin",
            email="superadmin@example.com",
            password="password123"
        )

        # Act
        created_admin = crud_user.create_superuser(db_session, user=admin_data)

        # Assert - Superusers should always be active (line 53 in crud_user.py)
        assert created_admin.is_active is True

    def test_create_superuser_hashes_password(self, db_session: Session):
        """Test that superuser password is also hashed."""
        # Arrange
        plain_password = "admin_secret_password"
        admin_data = UserCreate(
            username="secureadmin",
            email="secureadmin@example.com",
            password=plain_password
        )

        # Act
        created_admin = crud_user.create_superuser(db_session, user=admin_data)

        # Assert
        assert created_admin.hashed_password != plain_password
        assert created_admin.hashed_password.startswith("$2b$")
        assert verify_password(plain_password, created_admin.hashed_password) is True


@pytest.mark.unit
@pytest.mark.crud
class TestReadOperations:
    """Test read operations: get_user, get_user_by_email, get_users_count."""

    def test_get_user_by_id_success(self, db_session: Session, sample_user: models.User):
        """Test fetching existing user by ID."""
        # Act
        fetched_user = crud_user.get_user(db_session, id=sample_user.id)

        # Assert
        assert fetched_user is not None
        assert fetched_user.id == sample_user.id
        assert fetched_user.email == sample_user.email
        assert fetched_user.username == sample_user.username

    def test_get_user_by_id_not_found(self, db_session: Session):
        """Test that get_user returns None for nonexistent user ID."""
        # Arrange
        nonexistent_id = 99999

        # Act
        fetched_user = crud_user.get_user(db_session, id=nonexistent_id)

        # Assert
        assert fetched_user is None

    def test_get_user_by_email_success(self, db_session: Session, sample_user: models.User):
        """Test fetching user by email address (used for login)."""
        # Act
        fetched_user = crud_user.get_user_by_email(db_session, email=sample_user.email)

        # Assert
        assert fetched_user is not None
        assert fetched_user.id == sample_user.id
        assert fetched_user.email == sample_user.email

    def test_get_user_by_email_not_found(self, db_session: Session):
        """Test that get_user_by_email returns None for nonexistent email."""
        # Arrange
        nonexistent_email = "doesnotexist@example.com"

        # Act
        fetched_user = crud_user.get_user_by_email(db_session, email=nonexistent_email)

        # Assert
        assert fetched_user is None

    def test_get_users_count_empty_db(self, db_session: Session):
        """Test user count on empty database returns 0."""
        # Act
        count = crud_user.get_users_count(db_session)

        # Assert
        assert count == 0

    def test_get_users_count_multiple_users(self, db_session: Session, user_factory):
        """Test accurate user count with multiple users."""
        # Arrange - Create 5 users
        for i in range(5):
            user_factory(
                username=f"user{i}",
                email=f"user{i}@example.com"
            )

        # Act
        count = crud_user.get_users_count(db_session)

        # Assert
        assert count == 5


@pytest.mark.unit
@pytest.mark.crud
class TestUpdateUser:
    """Test update_user() function for user modifications.

    Note: Currently only supports password updates (line 20-26 in crud_user.py).
    """

    def test_update_user_password(self, db_session: Session, sample_user: models.User):
        """Test successful password update with new hash."""
        # Arrange
        old_hashed_password = sample_user.hashed_password
        new_password = "new_secure_password_456"
        update_data = UserUpdate(password=new_password)

        # Act
        updated_user = crud_user.update_user(
            db_session,
            user_id=sample_user.id,
            user=update_data
        )

        # Assert
        assert updated_user.id == sample_user.id
        assert updated_user.hashed_password != old_hashed_password
        assert verify_password(new_password, updated_user.hashed_password) is True

    def test_update_user_password_verification(self, db_session: Session, sample_user: models.User):
        """Test that old password no longer works after update."""
        # Arrange
        old_password = "testpassword123"  # From conftest.py sample_user
        new_password = "completely_new_password"
        update_data = UserUpdate(password=new_password)

        # Act
        updated_user = crud_user.update_user(
            db_session,
            user_id=sample_user.id,
            user=update_data
        )

        # Assert
        assert verify_password(old_password, updated_user.hashed_password) is False
        assert verify_password(new_password, updated_user.hashed_password) is True

    def test_update_user_returns_updated_model(self, db_session: Session, sample_user: models.User):
        """Test that update_user returns refreshed User object."""
        # Arrange
        update_data = UserUpdate(password="newpassword789")

        # Act
        updated_user = crud_user.update_user(
            db_session,
            user_id=sample_user.id,
            user=update_data
        )

        # Assert - Should return User model, not None
        assert isinstance(updated_user, models.User)
        assert updated_user.id == sample_user.id
        # Verify it's the same object reference from DB
        assert updated_user.email == sample_user.email


@pytest.mark.unit
@pytest.mark.crud
class TestAuthentication:
    """Test authentication and authorization helper functions.

    Functions tested: authenticate(), is_active(), is_superuser()
    """

    def test_authenticate_success(self, db_session: Session, sample_user: models.User):
        """Test successful authentication with valid credentials."""
        # Arrange
        email = sample_user.email
        password = "testpassword123"  # From conftest.py

        # Act
        authenticated_user = crud_user.authenticate(
            db_session,
            email=email,
            password=password
        )

        # Assert
        assert authenticated_user is not None
        assert authenticated_user.id == sample_user.id
        assert authenticated_user.email == email

    def test_authenticate_invalid_email(self, db_session: Session):
        """Test authenticate returns None for nonexistent email."""
        # Arrange
        nonexistent_email = "fake@example.com"
        password = "anypassword"

        # Act
        authenticated_user = crud_user.authenticate(
            db_session,
            email=nonexistent_email,
            password=password
        )

        # Assert
        assert authenticated_user is None

    def test_authenticate_invalid_password(self, db_session: Session, sample_user: models.User):
        """Test authenticate returns None for incorrect password."""
        # Arrange
        email = sample_user.email
        wrong_password = "wrongpassword"

        # Act
        authenticated_user = crud_user.authenticate(
            db_session,
            email=email,
            password=wrong_password
        )

        # Assert
        assert authenticated_user is None

    def test_is_active_true(self, db_session: Session, sample_user: models.User):
        """Test is_active returns True for active user."""
        # Act
        result = crud_user.is_active(sample_user)

        # Assert
        assert result is True

    def test_is_superuser_true(self, db_session: Session, sample_admin_user: models.User):
        """Test is_superuser returns True for admin user."""
        # Act
        result = crud_user.is_superuser(sample_admin_user)

        # Assert
        assert result is True
