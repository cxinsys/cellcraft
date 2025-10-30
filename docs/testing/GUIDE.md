# CellCraft Backend Testing - Test Writing Guide

> Comprehensive guide to writing high-quality unit tests for CellCraft backend

**Last Updated**: 2025-01-27
**Target Audience**: Backend developers, QA engineers
**Difficulty**: Intermediate

---

## Table of Contents

1. [Testing Philosophy](#testing-philosophy)
2. [AAA Pattern](#aaa-pattern)
3. [Test Naming Conventions](#test-naming-conventions)
4. [Fixture Usage](#fixture-usage)
5. [Mocking Strategies](#mocking-strategies)
6. [Pytest Markers](#pytest-markers)
7. [Assertion Best Practices](#assertion-best-practices)
8. [Coverage Targets](#coverage-targets)
9. [Common Pitfalls](#common-pitfalls)
10. [Code Examples](#code-examples)

---

## Testing Philosophy

### Core Principles

**1. Tests as Documentation**
- Tests should clearly demonstrate how code is intended to be used
- Test names should be self-documenting
- Comments should explain "why", not "what"

**2. Test Independence**
- Each test should be runnable in isolation
- No dependencies between tests
- Clean state before and after each test

**3. Comprehensive Coverage**
- **Happy Path**: Normal, expected usage
- **Edge Cases**: Boundary conditions, unusual inputs
- **Error Cases**: Invalid inputs, failure scenarios

**4. Maintainability**
- Keep tests simple and readable
- Avoid complex logic in tests
- DRY principle: Use fixtures and helpers

---

## AAA Pattern

### Overview

The **Arrange-Act-Assert (AAA)** pattern is the standard structure for writing clear, maintainable tests.

```python
def test_example():
    # Arrange: Set up test data and preconditions
    user_data = UserCreate(username="test", email="test@example.com", password="pass123")

    # Act: Execute the function under test
    result = crud_user.create_user(db_session, user=user_data)

    # Assert: Verify the expected outcome
    assert result.id is not None
    assert result.username == "test"
    assert result.email == "test@example.com"
```

### Detailed Breakdown

#### Arrange Phase

**Purpose**: Set up all preconditions and test data

**What to include**:
- Input data creation
- Mock setup
- Database state preparation
- Fixture initialization

**Example**:
```python
# Arrange
user_data = UserCreate(
    username="newuser",
    email="newuser@example.com",
    password="securepassword123"
)
# Setup complete - all inputs ready
```

#### Act Phase

**Purpose**: Execute the single operation being tested

**Guidelines**:
- Only ONE action per test
- Should be a single line (or minimal lines)
- No complex logic

**Example**:
```python
# Act
created_user = crud_user.create_user(db_session, user=user_data)
```

#### Assert Phase

**Purpose**: Verify the outcome matches expectations

**Guidelines**:
- Multiple assertions are OK if they verify related outcomes
- Use specific assertions (not just `assert result`)
- Include negative assertions when relevant

**Example**:
```python
# Assert
assert created_user.id is not None
assert created_user.username == "newuser"
assert created_user.email == "newuser@example.com"
assert created_user.hashed_password is not None
assert created_user.is_active is True
assert created_user.is_superuser is False
```

---

## Test Naming Conventions

### Structure

```
test_<function_name>_<scenario>_<expected_result>
```

### Examples

✅ **Good Names**:
```python
def test_create_user_success()
def test_create_user_with_duplicate_email_raises_error()
def test_authenticate_with_invalid_password_returns_none()
def test_get_user_by_id_not_found()
```

❌ **Bad Names**:
```python
def test_create()  # Too vague
def test_user()  # What about the user?
def test_1()  # No context
def test_error()  # What kind of error?
```

### Naming Guidelines

1. **Be Specific**: Name should describe exactly what is being tested
2. **Include Expected Outcome**: Make the expected result clear
3. **Use Underscores**: Separate words for readability
4. **Avoid Abbreviations**: Unless universally understood (e.g., "id", "api")
5. **Length is OK**: Descriptive names are better than short, cryptic ones

---

## Fixture Usage

### Available Fixtures (conftest.py)

#### db_session

**Purpose**: Provides isolated database session with automatic rollback

**Usage**:
```python
def test_create_user_success(db_session: Session):
    user_data = UserCreate(username="test", email="test@example.com", password="pass")
    created_user = crud_user.create_user(db_session, user=user_data)
    assert created_user.id is not None
```

**Behavior**:
- Creates new session for each test
- Automatically rolls back changes after test
- Ensures test isolation

---

#### sample_user

**Purpose**: Provides a pre-created test user

**Usage**:
```python
def test_get_user_by_id(db_session: Session, sample_user: models.User):
    fetched_user = crud_user.get_user(db_session, id=sample_user.id)
    assert fetched_user.id == sample_user.id
```

**Properties**:
```python
username = "testuser"
email = "testuser@example.com"
password = "testpassword123" (hashed)
is_active = True
is_superuser = False
```

---

#### sample_admin_user

**Purpose**: Provides a pre-created admin user

**Usage**:
```python
def test_is_superuser_true(sample_admin_user: models.User):
    result = crud_user.is_superuser(sample_admin_user)
    assert result is True
```

**Properties**:
```python
username = "adminuser"
email = "admin@example.com"
is_superuser = True
is_active = True
```

---

#### sample_plugin

**Purpose**: Provides a pre-created plugin for testing

**Usage**:
```python
def test_user_plugin_association(db_session: Session, sample_plugin: models.Plugin):
    user = create_user(db_session, user_data)
    db_session.refresh(user)
    assert len(user.plugins) == 1
    assert user.plugins[0].id == sample_plugin.id
```

---

#### user_factory

**Purpose**: Factory function to create multiple users dynamically

**Usage**:
```python
def test_create_multiple_users(db_session: Session, user_factory):
    # Create 5 users
    for i in range(5):
        user_factory(
            username=f"user{i}",
            email=f"user{i}@example.com"
        )

    count = crud_user.get_users_count(db_session)
    assert count == 5
```

---

### Custom Fixtures

**Creating Test-Specific Fixtures**:

```python
import pytest

@pytest.fixture
def sample_workflow(db_session: Session, sample_user: models.User) -> models.Workflow:
    """Create a test workflow for testing."""
    workflow = models.Workflow(
        name="Test Workflow",
        user_id=sample_user.id,
        json_data={"nodes": [], "connections": []}
    )
    db_session.add(workflow)
    db_session.commit()
    db_session.refresh(workflow)
    return workflow
```

---

## Mocking Strategies

### When to Mock

✅ **Mock These**:
- External API calls (HTTP requests)
- File system operations (when not testing I/O)
- Time-dependent functions (datetime.now)
- Expensive operations (large computations)
- Third-party services (email, payment gateways)

❌ **Don't Mock These**:
- Database operations (use test database)
- Your own business logic
- Simple utility functions
- Models and schemas

---

### Mocking Examples

#### Mock External API Call

```python
from unittest.mock import patch, MagicMock

def test_fetch_user_data_from_api():
    # Arrange
    mock_response = MagicMock()
    mock_response.json.return_value = {"id": 1, "name": "Test User"}
    mock_response.status_code = 200

    # Act
    with patch('requests.get', return_value=mock_response):
        result = fetch_user_data(user_id=1)

    # Assert
    assert result["name"] == "Test User"
```

#### Mock Datetime

```python
from datetime import datetime
from unittest.mock import patch

def test_create_timestamp():
    # Arrange
    fake_now = datetime(2025, 1, 27, 12, 0, 0)

    # Act
    with patch('datetime.datetime') as mock_datetime:
        mock_datetime.now.return_value = fake_now
        result = create_record()

    # Assert
    assert result.created_at == fake_now
```

#### Mock Environment Variables

```python
def test_config_with_custom_env(monkeypatch):
    # Arrange
    monkeypatch.setenv("POSTGRES_PORT", "5434")

    # Act
    settings = Settings()

    # Assert
    assert settings.POSTGRES_PORT == "5434"
```

---

## Pytest Markers

### Purpose

Markers allow you to categorize and selectively run tests.

### Available Markers

#### @pytest.mark.unit

**Usage**: Mark unit tests (isolated, fast tests)

```python
@pytest.mark.unit
def test_password_hashing():
    hashed = get_password_hash("plaintext")
    assert verify_password("plaintext", hashed) is True
```

**Run only unit tests**:
```bash
pytest tests/unit/ -m unit -v
```

---

#### @pytest.mark.integration

**Usage**: Mark integration tests (multiple components)

```python
@pytest.mark.integration
def test_full_user_registration_flow(db_session):
    # Register → Login → Fetch profile
    ...
```

**Run only integration tests**:
```bash
pytest tests/ -m integration -v
```

---

#### @pytest.mark.auth

**Usage**: Mark authentication-related tests

```python
@pytest.mark.auth
def test_login_with_valid_credentials():
    ...
```

---

#### @pytest.mark.crud

**Usage**: Mark CRUD operation tests

```python
@pytest.mark.crud
def test_create_user_success():
    ...
```

---

#### @pytest.mark.api

**Usage**: Mark API endpoint tests

```python
@pytest.mark.api
def test_register_endpoint_returns_201():
    ...
```

---

### Combining Markers

```python
@pytest.mark.unit
@pytest.mark.auth
class TestAuthentication:
    def test_authenticate_success(self):
        ...
```

**Run specific combination**:
```bash
pytest tests/ -m "unit and auth" -v
```

---

## Assertion Best Practices

### Use Specific Assertions

✅ **Good**:
```python
assert user.id is not None
assert user.username == "testuser"
assert user.is_active is True
```

❌ **Bad**:
```python
assert user  # Too vague
assert user.username  # What about the username?
```

---

### Assert Multiple Related Properties

```python
# Good: Verify all important properties
assert created_user.id is not None
assert created_user.username == "newuser"
assert created_user.email == "newuser@example.com"
assert created_user.hashed_password is not None
assert created_user.is_active is True
```

---

### Use pytest.raises for Exceptions

```python
from pydantic import ValidationError

def test_user_create_missing_email_raises_error():
    # Arrange
    invalid_data = {"username": "test", "password": "pass"}

    # Act & Assert
    with pytest.raises(ValidationError) as exc_info:
        UserCreate(**invalid_data)

    # Verify error details
    errors = exc_info.value.errors()
    assert any(error["loc"] == ("email",) for error in errors)
```

---

### Negative Assertions

**Verify what should NOT happen**:

```python
def test_password_not_stored_in_plaintext(db_session):
    # Arrange
    plain_password = "my_secret_password"
    user_data = UserCreate(username="user", email="user@example.com", password=plain_password)

    # Act
    created_user = crud_user.create_user(db_session, user=user_data)

    # Assert - Password should NOT be stored in plaintext
    assert created_user.hashed_password != plain_password
    assert created_user.hashed_password.startswith("$2b$")  # bcrypt format
```

---

## Coverage Targets

### Overall Targets

| Phase | Target Coverage | Current |
|-------|----------------|---------|
| Phase 1 | 20% | ✅ 22% |
| Phase 2 | 30% | 🔄 22% |
| Phase 3 | 50% | 📅 Planned |
| Phase 4 | 75% | 📅 Planned |

### Module-Level Targets

**Critical Modules** (90%+ coverage):
- `app/common/security.py` ✅ 100%
- `app/database/crud/crud_user.py` ✅ 100%
- `app/common/config.py` ✅ 93%

**Important Modules** (70%+ coverage):
- `app/database/models.py` 🔄 78%
- `app/common/error_utils.py` ✅ 89%
- `app/common/datatable_utils.py` ✅ 93%

**Secondary Modules** (50%+ coverage):
- Route modules (API endpoints)
- Utility modules

---

### Measuring Coverage

```bash
# Generate coverage report
pytest tests/unit/ --cov=app --cov-report=html

# View detailed report
open htmlcov/index.html
```

---

## Common Pitfalls

### ❌ Pitfall 1: Testing Implementation, Not Behavior

**Bad**:
```python
def test_create_user_calls_hash_function():
    # Testing internal implementation detail
    with patch('app.common.security.get_password_hash') as mock_hash:
        create_user(db_session, user_data)
        assert mock_hash.called
```

**Good**:
```python
def test_create_user_stores_hashed_password():
    # Testing observable behavior
    plain_password = "test123"
    user_data = UserCreate(username="user", email="user@example.com", password=plain_password)
    created_user = create_user(db_session, user_data)

    assert created_user.hashed_password != plain_password
    assert verify_password(plain_password, created_user.hashed_password) is True
```

---

### ❌ Pitfall 2: Over-Mocking

**Bad**:
```python
def test_create_user_with_too_many_mocks():
    # Mocking everything makes test meaningless
    with patch('crud_user.create_user') as mock_create:
        mock_create.return_value = MagicMock(id=1)
        result = crud_user.create_user(db_session, user_data)
        assert result.id == 1  # Testing the mock, not the code!
```

**Good**:
```python
def test_create_user_success(db_session):
    # Use real database operations
    user_data = UserCreate(username="user", email="user@example.com", password="pass")
    created_user = crud_user.create_user(db_session, user=user_data)
    assert created_user.id is not None
```

---

### ❌ Pitfall 3: Test Dependencies

**Bad**:
```python
# test_user.py
created_user_id = None

def test_1_create_user():
    global created_user_id
    user = create_user(...)
    created_user_id = user.id  # ❌ Shared state

def test_2_update_user():
    # ❌ Depends on test_1
    user = get_user(created_user_id)
    ...
```

**Good**:
```python
def test_create_user(db_session):
    # Independent test
    user = create_user(db_session, ...)
    assert user.id is not None

def test_update_user(db_session, sample_user):
    # Independent test with fixture
    updated = update_user(db_session, sample_user.id, ...)
    assert updated.username == "newname"
```

---

### ❌ Pitfall 4: Unclear Test Names

**Bad**:
```python
def test_user():  # What about the user?
def test_error():  # What kind of error?
def test_success():  # Success of what?
```

**Good**:
```python
def test_create_user_with_duplicate_email_raises_error()
def test_authenticate_with_invalid_password_returns_none()
def test_update_user_password_success()
```

---

### ❌ Pitfall 5: Testing Multiple Things

**Bad**:
```python
def test_user_operations():
    # ❌ Testing create, update, delete in one test
    user = create_user(...)
    assert user.id is not None

    updated = update_user(...)
    assert updated.username == "new"

    delete_user(...)
    assert get_user(...) is None
```

**Good**:
```python
def test_create_user_success():
    user = create_user(...)
    assert user.id is not None

def test_update_user_success():
    user = update_user(...)
    assert user.username == "new"

def test_delete_user_success():
    delete_user(...)
    assert get_user(...) is None
```

---

## Code Examples

### Example 1: Simple Function Test

```python
from app.common.security import verify_password, get_password_hash

@pytest.mark.unit
def test_password_verification_success():
    """Test that password verification works with correct password."""
    # Arrange
    plain_password = "my_secure_password"
    hashed_password = get_password_hash(plain_password)

    # Act
    result = verify_password(plain_password, hashed_password)

    # Assert
    assert result is True
```

---

### Example 2: Database CRUD Test

```python
from app.database.crud import crud_user
from app.database.schemas.user import UserCreate

@pytest.mark.unit
@pytest.mark.crud
def test_create_user_success(db_session: Session):
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
```

---

### Example 3: Schema Validation Test

```python
from pydantic import ValidationError
from app.database.schemas.user import UserCreate

@pytest.mark.unit
@pytest.mark.api
def test_user_create_missing_email_raises_validation_error():
    """Test UserCreate rejects missing required email field."""
    # Arrange
    invalid_data = {
        "username": "testuser",
        "password": "securepassword123"
        # Missing required email field
    }

    # Act & Assert
    with pytest.raises(ValidationError) as exc_info:
        UserCreate(**invalid_data)

    # Verify error details
    errors = exc_info.value.errors()
    assert any(error["loc"] == ("email",) for error in errors)
    assert any("field required" in error["msg"].lower() or
               "missing" in error["msg"].lower()
               for error in errors)
```

---

### Example 4: Testing Edge Cases

```python
@pytest.mark.unit
@pytest.mark.crud
def test_create_user_with_zero_plugins(db_session: Session):
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
```

---

### Example 5: Testing Error Cases

```python
@pytest.mark.unit
@pytest.mark.crud
def test_get_user_by_id_not_found(db_session: Session):
    """Test that get_user returns None for nonexistent user ID."""
    # Arrange
    nonexistent_id = 99999

    # Act
    fetched_user = crud_user.get_user(db_session, id=nonexistent_id)

    # Assert
    assert fetched_user is None
```

---

### Example 6: Using Fixtures

```python
@pytest.mark.unit
@pytest.mark.crud
def test_get_user_by_email_success(db_session: Session, sample_user: models.User):
    """Test fetching user by email address (used for login)."""
    # Act
    fetched_user = crud_user.get_user_by_email(db_session, email=sample_user.email)

    # Assert
    assert fetched_user is not None
    assert fetched_user.id == sample_user.id
    assert fetched_user.email == sample_user.email
```

---

### Example 7: Testing Complex Logic

```python
@pytest.mark.unit
def test_assemble_cors_origins_with_json_array_string():
    """Test CORS origins validator with JSON array string - verify actual parsing."""
    # Arrange
    cors_json = '["http://localhost:8080", "http://localhost:3000"]'

    # Act
    result = Settings.assemble_cors_origins(cors_json)

    # Assert - Result should be a list, not a string
    assert isinstance(result, list), "JSON string should be parsed to list"
    assert result == ["http://localhost:8080", "http://localhost:3000"]
    assert len(result) == 2
```

---

## Next Steps

1. **Read Codebase Analysis**: See [CODEBASE-ANALYSIS.md](./CODEBASE-ANALYSIS.md) for complexity insights
2. **Review Progress**: See [PROGRESS.md](./PROGRESS.md) for current status
3. **Set Up Environment**: See [SETUP.md](./SETUP.md) if not already configured
4. **Start Writing Tests**: Apply these principles to new test development

---

## Quick Reference

### Test Structure Checklist

- [ ] Clear, descriptive test name
- [ ] AAA pattern (Arrange, Act, Assert)
- [ ] Proper pytest markers
- [ ] Independent from other tests
- [ ] Covers happy path, edge cases, and errors
- [ ] Uses appropriate fixtures
- [ ] Assertions are specific and comprehensive
- [ ] No complex logic in test
- [ ] Follows project conventions

---

**Happy Testing! 🧪**
