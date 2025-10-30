"""
Unit tests for configuration settings.

Test coverage:
- Celery task routing
- CORS origins validation
- Settings initialization
- Configuration caching
"""
import pytest
import json

from app.common.config import route_task, Settings, get_settings


@pytest.mark.unit
@pytest.mark.config
class TestCeleryTaskRouting:
    """Test Celery task routing function."""

    def test_route_task_with_queue_prefix(self):
        """Test task routing when task name contains queue prefix."""
        # Arrange
        task_name = "workflow_task:process_data"

        # Act
        result = route_task(task_name, args=[], kwargs={}, options={})

        # Assert
        assert result == {"queue": "workflow_task"}

    def test_route_task_without_queue_prefix(self):
        """Test task routing when task name has no queue prefix."""
        # Arrange
        task_name = "simple_task"

        # Act
        result = route_task(task_name, args=[], kwargs={}, options={})

        # Assert
        assert result == {"queue": "celery"}


@pytest.mark.unit
@pytest.mark.config
class TestCorsOriginsValidator:
    """Test CORS origins validation."""

    def test_assemble_cors_origins_from_comma_separated_string(self):
        """Test CORS origins validator with comma-separated string."""
        # Arrange
        cors_string = "http://localhost:8080, http://localhost:3000"

        # Act
        result = Settings.assemble_cors_origins(cors_string)

        # Assert
        assert result == ["http://localhost:8080", "http://localhost:3000"]

    def test_assemble_cors_origins_from_list(self):
        """Test CORS origins validator with list input."""
        # Arrange
        cors_list = ["http://localhost:8080", "http://localhost:3000"]

        # Act
        result = Settings.assemble_cors_origins(cors_list)

        # Assert
        assert result == cors_list

    def test_assemble_cors_origins_with_json_array_string(self):
        """Test CORS origins validator with JSON array string - verify actual parsing."""
        # Arrange
        cors_json = '["http://localhost:8080", "http://localhost:3000"]'

        # Act
        result = Settings.assemble_cors_origins(cors_json)

        # Assert - Result should be a list, not a string
        assert isinstance(result, list), "JSON string should be parsed to list"
        assert result == ["http://localhost:8080", "http://localhost:3000"]
        assert len(result) == 2

    def test_assemble_cors_origins_with_invalid_json_raises_error(self):
        """Test CORS origins validator with invalid JSON raises appropriate error."""
        # Arrange
        invalid_json = '["http://localhost:8080", invalid]'  # Invalid JSON

        # Act & Assert
        with pytest.raises((json.JSONDecodeError, ValueError)):
            Settings.assemble_cors_origins(invalid_json)


@pytest.mark.unit
@pytest.mark.config
class TestSettingsInitialization:
    """Test Settings class initialization."""

    def test_settings_database_uri_construction(self):
        """Test database URI is correctly constructed from class-level environment variables."""
        # Act - Settings class reads environment at class definition time
        settings = Settings()

        # Assert - Verify URI components are not None (loaded from environment)
        assert settings.POSTGRES_USER is not None, "POSTGRES_USER should be loaded from environment"
        assert settings.POSTGRES_PASSWORD is not None, "POSTGRES_PASSWORD should be loaded from environment"
        assert settings.POSTGRES_HOST is not None, "POSTGRES_HOST should be loaded from environment"
        assert settings.POSTGRES_PORT is not None, "POSTGRES_PORT should be loaded from environment"
        assert settings.POSTGRES_DB is not None, "POSTGRES_DB should be loaded from environment"

        # Verify URI is correctly constructed from loaded values
        expected_uri = f"postgresql://{settings.POSTGRES_USER}:{settings.POSTGRES_PASSWORD}@{settings.POSTGRES_HOST}:{settings.POSTGRES_PORT}/{settings.POSTGRES_DB}"
        assert settings.SQLALCHEMY_DATABASE_URI == expected_uri

        # Verify URI format and structure
        assert settings.SQLALCHEMY_DATABASE_URI.startswith("postgresql://")
        assert settings.POSTGRES_DB in settings.SQLALCHEMY_DATABASE_URI
        assert settings.POSTGRES_USER in settings.SQLALCHEMY_DATABASE_URI
        assert settings.POSTGRES_HOST in settings.SQLALCHEMY_DATABASE_URI

    def test_settings_database_uri_construction_with_monkeypatch(self, monkeypatch):
        """Test database URI construction with isolated environment variables."""
        # Arrange - Set known environment variables
        monkeypatch.setenv("POSTGRES_USER", "test_user")
        monkeypatch.setenv("POSTGRES_PASSWORD", "test_pass")
        monkeypatch.setenv("POSTGRES_HOST", "localhost")
        monkeypatch.setenv("POSTGRES_PORT", "5432")
        monkeypatch.setenv("POSTGRES_DB", "test_db")

        # Act
        settings = Settings()

        # Assert
        expected_uri = "postgresql://test_user:test_pass@localhost:5432/test_db"
        assert settings.SQLALCHEMY_DATABASE_URI == expected_uri


@pytest.mark.unit
@pytest.mark.config
class TestGetSettings:
    """Test get_settings function and caching behavior."""

    def test_get_settings_returns_development_config(self):
        """Test get_settings returns DevelopmentConfig by default."""
        # Act
        settings = get_settings()

        # Assert
        assert settings is not None
        assert hasattr(settings, 'SQLALCHEMY_DATABASE_URI')
        assert hasattr(settings, 'CELERY_BROKER_URL')
        assert settings.PROJECT_NAME == "test"

    def test_get_settings_returns_same_instance(self):
        """Test get_settings caching returns the same instance."""
        # Act
        settings1 = get_settings()
        settings2 = get_settings()

        # Assert - Verify it's the same object instance
        assert settings1 is settings2
        assert id(settings1) == id(settings2)

    def test_get_settings_caching_preserves_attributes(self):
        """Test get_settings caching preserves attributes across calls."""
        # Arrange & Act
        settings1 = get_settings()
        project_name = settings1.PROJECT_NAME
        database_uri = settings1.SQLALCHEMY_DATABASE_URI

        settings2 = get_settings()

        # Assert - Verify cached instance has same attributes
        assert settings2.PROJECT_NAME == project_name
        assert settings2.SQLALCHEMY_DATABASE_URI == database_uri
