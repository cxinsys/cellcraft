"""
Unit tests for SQLAlchemy database models.

Test coverage:
- User model: field defaults, constraints, relationships
- Plugin model: JSONB fields, enum validation, timestamps, many-to-many
- Workflow model: JSONB structure, foreign keys
- File model: constraints, relationships
- Task model: multiple foreign keys
- Cross-model relationships and cascade behaviors

Phase 2.2: SQLAlchemy Model Testing
Quality Target: 8/10 with AAA pattern compliance
"""
import pytest
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from app.database import models
from app.common.enums import PluginType


@pytest.mark.unit
@pytest.mark.models
class TestUserModel:
    """Test User model field validation, constraints, and relationships."""

    def test_user_model_field_defaults(self, db_session: Session):
        """Test User model field default values are set correctly."""
        # Arrange
        user = models.User(
            username="testuser",
            email="test@example.com",
            hashed_password="hashed_password_123"
        )

        # Act
        db_session.add(user)
        db_session.commit()
        db_session.refresh(user)

        # Assert - Verify default values
        assert user.is_active is True, "is_active should default to True"
        assert user.is_superuser is False, "is_superuser should default to False"
        assert user.id is not None, "id should be auto-generated"

    def test_user_email_unique_constraint(self, db_session: Session):
        """Test User email unique constraint enforcement."""
        # Arrange
        user1 = models.User(
            username="user1",
            email="duplicate@example.com",
            hashed_password="hash1"
        )
        db_session.add(user1)
        db_session.commit()

        # Act & Assert - Attempt to create duplicate email
        user2 = models.User(
            username="user2",
            email="duplicate@example.com",  # Duplicate email
            hashed_password="hash2"
        )
        db_session.add(user2)

        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()

        assert "unique constraint" in str(exc_info.value).lower() or "duplicate key" in str(exc_info.value).lower()
        db_session.rollback()  # Clean up failed transaction

    def test_user_nullable_constraints(self, db_session: Session):
        """Test User model nullable constraints for required fields."""
        # Test 1: username cannot be NULL
        # Arrange
        user_no_username = models.User(
            username=None,  # NULL username
            email="test@example.com",
            hashed_password="hash"
        )
        db_session.add(user_no_username)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 2: email cannot be NULL
        # Arrange
        user_no_email = models.User(
            username="testuser",
            email=None,  # NULL email
            hashed_password="hash"
        )
        db_session.add(user_no_email)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 3: hashed_password cannot be NULL
        # Arrange
        user_no_password = models.User(
            username="testuser",
            email="test@example.com",
            hashed_password=None  # NULL password
        )
        db_session.add(user_no_password)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

    def test_user_hashed_password_not_plaintext(self, db_session: Session, sample_user: models.User):
        """Test User hashed_password is not stored as plaintext."""
        # Arrange
        plaintext_password = "testpassword123"  # Known plaintext from sample_user fixture

        # Act
        db_session.refresh(sample_user)

        # Assert
        assert sample_user.hashed_password != plaintext_password, \
            "Password should be hashed, not stored as plaintext"
        assert len(sample_user.hashed_password) > len(plaintext_password), \
            "Hashed password should be longer than plaintext"
        assert sample_user.hashed_password.startswith("$2b$"), \
            "bcrypt hash should start with $2b$ algorithm identifier"

    def test_user_relationships_files(self, db_session: Session, sample_user: models.User):
        """Test User.files relationship loads correctly (one-to-many)."""
        # Arrange
        file1 = models.File(
            file_name="test1.h5ad",
            file_size="1024",
            file_path="/app/user_data/testuser/test1.h5ad",
            folder="/app/user_data/testuser",
            user_id=sample_user.id
        )
        file2 = models.File(
            file_name="test2.h5ad",
            file_size="2048",
            file_path="/app/user_data/testuser/test2.h5ad",
            folder="/app/user_data/testuser",
            user_id=sample_user.id
        )
        db_session.add_all([file1, file2])
        db_session.commit()

        # Act
        db_session.refresh(sample_user)

        # Assert - Verify relationship loads
        assert len(sample_user.files) == 2, "User should have 2 files"
        assert file1 in sample_user.files, "file1 should be in user.files"
        assert file2 in sample_user.files, "file2 should be in user.files"

        # Assert - Verify bidirectional relationship
        assert file1.user == sample_user, "file1.user should reference sample_user"
        assert file2.user == sample_user, "file2.user should reference sample_user"

    def test_user_relationships_workflows(self, db_session: Session, sample_user: models.User):
        """Test User.workflows relationship loads correctly (one-to-many)."""
        # Arrange
        workflow1 = models.Workflow(
            title="Test Workflow 1",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=sample_user.id
        )
        workflow2 = models.Workflow(
            title="Test Workflow 2",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=sample_user.id
        )
        db_session.add_all([workflow1, workflow2])
        db_session.commit()

        # Act
        db_session.refresh(sample_user)

        # Assert - Verify relationship loads
        assert len(sample_user.workflows) == 2, "User should have 2 workflows"
        assert workflow1 in sample_user.workflows, "workflow1 should be in user.workflows"
        assert workflow2 in sample_user.workflows, "workflow2 should be in user.workflows"

        # Assert - Verify bidirectional relationship
        assert workflow1.user == sample_user, "workflow1.user should reference sample_user"
        assert workflow2.user == sample_user, "workflow2.user should reference sample_user"

    def test_user_relationships_tasks(self, db_session: Session, sample_user: models.User, sample_workflow: models.Workflow):
        """Test User.tasks relationship loads correctly (one-to-many)."""
        # Arrange
        from datetime import datetime
        task1 = models.Task(
            task_id="task_123",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id
        )
        task2 = models.Task(
            task_id="task_456",
            start_time=datetime.now(),
            status="completed",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id
        )
        db_session.add_all([task1, task2])
        db_session.commit()

        # Act
        db_session.refresh(sample_user)

        # Assert - Verify relationship loads
        assert len(sample_user.tasks) == 2, "User should have 2 tasks"
        assert task1 in sample_user.tasks, "task1 should be in user.tasks"
        assert task2 in sample_user.tasks, "task2 should be in user.tasks"

        # Assert - Verify bidirectional relationship
        assert task1.user == sample_user, "task1.user should reference sample_user"
        assert task2.user == sample_user, "task2.user should reference sample_user"

    def test_user_relationships_plugins_many_to_many(self, db_session: Session, sample_user: models.User, sample_plugin: models.Plugin):
        """Test User.plugins many-to-many relationship via association table."""
        # Arrange - Associate plugin with user
        sample_user.plugins.append(sample_plugin)
        db_session.commit()

        # Act
        db_session.refresh(sample_user)
        db_session.refresh(sample_plugin)

        # Assert - Verify User -> Plugins relationship
        assert len(sample_user.plugins) == 1, "User should have 1 associated plugin"
        assert sample_plugin in sample_user.plugins, "sample_plugin should be in user.plugins"

        # Assert - Verify Plugin -> Users relationship (bidirectional)
        assert len(sample_plugin.users) == 1, "Plugin should have 1 associated user"
        assert sample_user in sample_plugin.users, "sample_user should be in plugin.users"

        # Assert - Verify association table entry
        result = db_session.execute(
            "SELECT user_id, plugin_id FROM user_plugin_association WHERE user_id = :user_id AND plugin_id = :plugin_id",
            {"user_id": sample_user.id, "plugin_id": sample_plugin.id}
        ).fetchone()
        assert result is not None, "Association table should have entry for user-plugin relationship"


@pytest.mark.unit
@pytest.mark.models
class TestPluginModel:
    """Test Plugin model field validation, JSONB fields, enum, and timestamps."""

    def test_plugin_model_field_defaults(self, db_session: Session):
        """Test Plugin model field default values are set correctly."""
        # Arrange
        plugin = models.Plugin(
            name="TestPlugin",
            description="Test plugin description",
            author="Test Author",
            plugin_path="/plugin/local/TestPlugin/",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "rule test: output: 'test.txt'"}
        )

        # Act
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - Verify default values
        assert plugin.use_gpu is False, "use_gpu should default to False"
        assert plugin.source == "local", "source should default to 'local'"
        assert plugin.is_editable is True, "is_editable should default to True"
        assert plugin.id is not None, "id should be auto-generated"
        assert plugin.created_at is not None, "created_at should be auto-set"
        assert plugin.updated_at is not None, "updated_at should be auto-set"

    def test_plugin_nullable_constraints(self, db_session: Session):
        """Test Plugin model nullable constraints for required fields."""
        # Test 1: name cannot be NULL
        # Arrange
        plugin_no_name = models.Plugin(
            name=None,  # NULL name
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin_no_name)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 2: description cannot be NULL
        # Arrange
        plugin_no_description = models.Plugin(
            name="Test",
            description=None,  # NULL description
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin_no_description)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 3: author cannot be NULL
        # Arrange
        plugin_no_author = models.Plugin(
            name="Test",
            description="Test",
            author=None,  # NULL author
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin_no_author)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

    def test_plugin_enum_validation(self, db_session: Session):
        """Test Plugin plugin_type enum accepts valid PluginType values."""
        # Test 1: PluginType.ANALYSIS is valid
        # Arrange
        plugin_analysis = models.Plugin(
            name="TestAnalysis",
            description="Test",
            author="Test",
            plugin_path="/path",
            plugin_type=PluginType.ANALYSIS,  # Valid enum
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )

        # Act
        db_session.add(plugin_analysis)
        db_session.commit()
        db_session.refresh(plugin_analysis)

        # Assert
        assert plugin_analysis.plugin_type == PluginType.ANALYSIS
        assert plugin_analysis.plugin_type.value == "analysis"

        # Test 2: PluginType.VISUALIZATION is valid
        # Arrange
        plugin_viz = models.Plugin(
            name="TestVisualization",
            description="Test",
            author="Test",
            plugin_path="/path/viz",
            plugin_type=PluginType.VISUALIZATION,  # Valid enum
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )

        # Act
        db_session.add(plugin_viz)
        db_session.commit()
        db_session.refresh(plugin_viz)

        # Assert
        assert plugin_viz.plugin_type == PluginType.VISUALIZATION
        assert plugin_viz.plugin_type.value == "visualization"

    def test_plugin_invalid_enum_raises_error(self, db_session: Session):
        """Test Plugin plugin_type rejects invalid enum values."""
        # Arrange - Create plugin with invalid enum value
        # Note: This tests database-level constraint, not Python-level

        # Create plugin with valid enum first
        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            plugin_type=PluginType.ANALYSIS,
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin)
        db_session.commit()

        # Act & Assert - Try to update with invalid string directly at DB level
        with pytest.raises(Exception):  # Could be IntegrityError or DataError depending on DB
            db_session.execute(
                "UPDATE plugins SET plugin_type = 'invalid_type' WHERE id = :id",
                {"id": plugin.id}
            )
            db_session.commit()
        db_session.rollback()

    def test_plugin_jsonb_dependencies_structure(self, db_session: Session):
        """Test Plugin dependencies JSONB field accepts valid structure."""
        # Arrange - Valid dependencies structure from sample_plugin fixture
        valid_dependencies = {
            "requirements.txt": "numpy==1.24.0\npandas==2.0.0\nscanpy==1.9.0"
        }

        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            dependencies=valid_dependencies,  # JSONB field
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )

        # Act
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - JSONB structure preserved
        assert plugin.dependencies == valid_dependencies
        assert "requirements.txt" in plugin.dependencies
        assert "numpy==1.24.0" in plugin.dependencies["requirements.txt"]

    def test_plugin_jsonb_drawflow_structure(self, db_session: Session):
        """Test Plugin drawflow JSONB field accepts valid structure."""
        # Arrange - Valid drawflow structure from sample_plugin fixture
        valid_drawflow = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "name": "test_node",
                            "data": {
                                "parameters": [
                                    {"name": "param1", "defaultValue": "test"}
                                ]
                            }
                        }
                    }
                }
            }
        }

        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow=valid_drawflow,  # JSONB field
            rules={"Snakefile": "test"}
        )

        # Act
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - JSONB structure preserved
        assert plugin.drawflow == valid_drawflow
        assert "drawflow" in plugin.drawflow
        assert "Home" in plugin.drawflow["drawflow"]
        assert "data" in plugin.drawflow["drawflow"]["Home"]
        assert "1" in plugin.drawflow["drawflow"]["Home"]["data"]

    def test_plugin_jsonb_rules_structure(self, db_session: Session):
        """Test Plugin rules JSONB field accepts valid structure."""
        # Arrange - Valid rules structure from sample_plugin fixture
        valid_rules = {
            "Snakefile": "rule test:\n    output: 'test.txt'\n    shell: 'echo test > {output}'"
        }

        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules=valid_rules  # JSONB field
        )

        # Act
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - JSONB structure preserved
        assert plugin.rules == valid_rules
        assert "Snakefile" in plugin.rules
        assert "rule test:" in plugin.rules["Snakefile"]

    def test_plugin_timestamps_auto_set(self, db_session: Session):
        """Test Plugin created_at and updated_at are auto-set on creation."""
        # Arrange
        from datetime import datetime, timezone

        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )

        # Act
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - Timestamps are auto-set
        assert plugin.created_at is not None, "created_at should be auto-set"
        assert plugin.updated_at is not None, "updated_at should be auto-set"

        # Assert - Timestamps have timezone info
        assert plugin.created_at.tzinfo is not None, "created_at should be timezone-aware"
        assert plugin.updated_at.tzinfo is not None, "updated_at should be timezone-aware"

        # Assert - Initially, created_at == updated_at (within 1 second)
        time_diff = abs((plugin.updated_at - plugin.created_at).total_seconds())
        assert time_diff < 1, "updated_at should equal created_at on initial creation"

    def test_plugin_updated_at_auto_update(self, db_session: Session):
        """Test Plugin updated_at is automatically updated on modification."""
        # Arrange - Create plugin
        import time
        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        original_created_at = plugin.created_at
        original_updated_at = plugin.updated_at

        # Wait to ensure timestamp difference
        time.sleep(1)

        # Act - Update plugin
        plugin.description = "Updated description"
        db_session.commit()
        db_session.refresh(plugin)

        # Assert - created_at unchanged
        assert plugin.created_at == original_created_at, \
            "created_at should not change on update"

        # Assert - updated_at changed
        assert plugin.updated_at > original_updated_at, \
            "updated_at should be automatically updated on modification"

        # Assert - updated_at is later than created_at
        assert plugin.updated_at > plugin.created_at, \
            "updated_at should be later than created_at after modification"

    def test_plugin_many_to_many_users(self, db_session: Session, sample_user: models.User):
        """Test Plugin.users many-to-many relationship via association table."""
        # Arrange - Create plugin and associate with user
        plugin = models.Plugin(
            name="TestPlugin",
            description="Test",
            author="Test",
            plugin_path="/path",
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "test"}
        )
        db_session.add(plugin)
        db_session.commit()

        # Act - Associate plugin with user
        plugin.users.append(sample_user)
        db_session.commit()

        db_session.refresh(plugin)
        db_session.refresh(sample_user)

        # Assert - Plugin -> Users relationship
        assert len(plugin.users) == 1, "Plugin should have 1 associated user"
        assert sample_user in plugin.users, "sample_user should be in plugin.users"

        # Assert - User -> Plugins relationship (bidirectional)
        assert plugin in sample_user.plugins, "plugin should be in sample_user.plugins"

        # Assert - Verify association table entry
        result = db_session.execute(
            "SELECT user_id, plugin_id FROM user_plugin_association WHERE user_id = :user_id AND plugin_id = :plugin_id",
            {"user_id": sample_user.id, "plugin_id": plugin.id}
        ).fetchone()
        assert result is not None, "Association table should have entry for user-plugin relationship"


@pytest.mark.unit
@pytest.mark.models
class TestWorkflowModel:
    """Test Workflow model field validation, JSONB structure, and relationships."""

    def test_workflow_nullable_constraints(self, db_session: Session):
        """Test Workflow model nullable constraints for required fields."""
        # Test 1: title cannot be NULL
        # Arrange
        workflow_no_title = models.Workflow(
            title=None,  # NULL title
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=1
        )
        db_session.add(workflow_no_title)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 2: workflow_info cannot be NULL
        # Arrange
        workflow_no_info = models.Workflow(
            title="Test Workflow",
            workflow_info=None,  # NULL workflow_info
            user_id=1
        )
        db_session.add(workflow_no_info)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

    def test_workflow_jsonb_structure(self, db_session: Session, sample_user: models.User):
        """Test Workflow workflow_info JSONB field accepts valid structure."""
        # Arrange - Valid workflow_info structure with drawflow
        valid_workflow_info = {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "data_input",
                            "data": {},
                            "inputs": {},
                            "outputs": {}
                        },
                        "2": {
                            "id": 2,
                            "name": "analysis_node",
                            "data": {"plugin": "TENET"},
                            "inputs": {},
                            "outputs": {}
                        }
                    }
                }
            }
        }

        workflow = models.Workflow(
            title="Complex Workflow",
            workflow_info=valid_workflow_info,  # JSONB field
            user_id=sample_user.id
        )

        # Act
        db_session.add(workflow)
        db_session.commit()
        db_session.refresh(workflow)

        # Assert - JSONB structure preserved
        assert workflow.workflow_info == valid_workflow_info
        assert "drawflow" in workflow.workflow_info
        assert "Home" in workflow.workflow_info["drawflow"]
        assert "data" in workflow.workflow_info["drawflow"]["Home"]
        assert "1" in workflow.workflow_info["drawflow"]["Home"]["data"]
        assert "2" in workflow.workflow_info["drawflow"]["Home"]["data"]

    def test_workflow_foreign_key_user(self, db_session: Session):
        """Test Workflow user_id foreign key constraint enforcement."""
        # Arrange - Create workflow with invalid user_id
        invalid_user_id = 99999  # Non-existent user

        workflow = models.Workflow(
            title="Test Workflow",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=invalid_user_id
        )
        db_session.add(workflow)

        # Act & Assert - Foreign key constraint should fail
        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()

        assert "foreign key" in str(exc_info.value).lower() or "violates" in str(exc_info.value).lower()
        db_session.rollback()

    def test_workflow_relationship_user(self, db_session: Session, sample_user: models.User):
        """Test Workflow.user relationship loads correctly (many-to-one)."""
        # Arrange
        workflow = models.Workflow(
            title="Test Workflow",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()

        # Act
        db_session.refresh(workflow)
        db_session.refresh(sample_user)

        # Assert - Workflow -> User relationship
        assert workflow.user is not None, "workflow.user should be loaded"
        assert workflow.user == sample_user, "workflow.user should reference sample_user"
        assert workflow.user.id == sample_user.id

        # Assert - User -> Workflows relationship (bidirectional)
        assert workflow in sample_user.workflows, "workflow should be in sample_user.workflows"

    def test_workflow_relationship_tasks(self, db_session: Session, sample_user: models.User):
        """Test Workflow.tasks relationship loads correctly (one-to-many)."""
        # Arrange - Create workflow first
        workflow = models.Workflow(
            title="Test Workflow",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()

        # Create tasks associated with workflow
        from datetime import datetime
        task1 = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=workflow.id
        )
        task2 = models.Task(
            task_id="task_002",
            start_time=datetime.now(),
            status="completed",
            user_id=sample_user.id,
            workflow_id=workflow.id
        )
        db_session.add_all([task1, task2])
        db_session.commit()

        # Act
        db_session.refresh(workflow)

        # Assert - Workflow -> Tasks relationship
        assert len(workflow.tasks) == 2, "Workflow should have 2 tasks"
        assert task1 in workflow.tasks, "task1 should be in workflow.tasks"
        assert task2 in workflow.tasks, "task2 should be in workflow.tasks"

        # Assert - Task -> Workflow relationship (bidirectional)
        assert task1.workflows == workflow, "task1.workflows should reference workflow"
        assert task2.workflows == workflow, "task2.workflows should reference workflow"


@pytest.mark.unit
@pytest.mark.models
class TestFileModel:
    """Test File model field validation and relationships."""

    def test_file_nullable_constraints(self, db_session: Session):
        """Test File model nullable constraints for all required fields."""
        # Test 1: file_name cannot be NULL
        # Arrange
        file_no_name = models.File(
            file_name=None,  # NULL file_name
            file_size="1024",
            file_path="/path/to/file",
            folder="/path/to",
            user_id=1
        )
        db_session.add(file_no_name)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 2: file_size cannot be NULL
        # Arrange
        file_no_size = models.File(
            file_name="test.h5ad",
            file_size=None,  # NULL file_size
            file_path="/path/to/file",
            folder="/path/to",
            user_id=1
        )
        db_session.add(file_no_size)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 3: file_path cannot be NULL
        # Arrange
        file_no_path = models.File(
            file_name="test.h5ad",
            file_size="1024",
            file_path=None,  # NULL file_path
            folder="/path/to",
            user_id=1
        )
        db_session.add(file_no_path)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 4: folder cannot be NULL
        # Arrange
        file_no_folder = models.File(
            file_name="test.h5ad",
            file_size="1024",
            file_path="/path/to/file",
            folder=None,  # NULL folder
            user_id=1
        )
        db_session.add(file_no_folder)

        # Act & Assert
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

    def test_file_foreign_key_user(self, db_session: Session):
        """Test File user_id foreign key constraint enforcement."""
        # Arrange - Create file with invalid user_id
        invalid_user_id = 99999  # Non-existent user

        file_record = models.File(
            file_name="test.h5ad",
            file_size="1024",
            file_path="/path/to/test.h5ad",
            folder="/path/to",
            user_id=invalid_user_id
        )
        db_session.add(file_record)

        # Act & Assert - Foreign key constraint should fail
        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()

        assert "foreign key" in str(exc_info.value).lower() or "violates" in str(exc_info.value).lower()
        db_session.rollback()

    def test_file_relationship_user(self, db_session: Session, sample_user: models.User):
        """Test File.user relationship loads correctly (many-to-one)."""
        # Arrange
        file_record = models.File(
            file_name="test.h5ad",
            file_size="1048576",  # 1MB as string
            file_path="/app/user_data/testuser/test.h5ad",
            folder="/app/user_data/testuser",
            user_id=sample_user.id
        )
        db_session.add(file_record)
        db_session.commit()

        # Act
        db_session.refresh(file_record)
        db_session.refresh(sample_user)

        # Assert - File -> User relationship
        assert file_record.user is not None, "file.user should be loaded"
        assert file_record.user == sample_user, "file.user should reference sample_user"
        assert file_record.user.id == sample_user.id

        # Assert - User -> Files relationship (bidirectional)
        assert file_record in sample_user.files, "file should be in sample_user.files"

    def test_file_size_as_string(self, db_session: Session, sample_user: models.User):
        """Test File file_size is stored as String type, not Integer."""
        # Arrange - Create file with large file size as string
        large_size_str = "1073741824"  # 1GB as string

        file_record = models.File(
            file_name="large_dataset.h5ad",
            file_size=large_size_str,  # String type
            file_path="/app/user_data/testuser/large_dataset.h5ad",
            folder="/app/user_data/testuser",
            user_id=sample_user.id
        )

        # Act
        db_session.add(file_record)
        db_session.commit()
        db_session.refresh(file_record)

        # Assert - file_size is stored as string
        assert isinstance(file_record.file_size, str), "file_size should be String type"
        assert file_record.file_size == large_size_str
        assert file_record.file_size == "1073741824"

        # Assert - Can perform string operations
        assert file_record.file_size.isdigit(), "file_size string should contain only digits"
        assert int(file_record.file_size) == 1073741824, "file_size string should be convertible to int"


# ==============================================================================
# Task Model Tests
# ==============================================================================

@pytest.mark.unit
@pytest.mark.models
class TestTaskModel:
    """Test Task model field validation, FK relationships, and optional fields."""

    def test_task_nullable_constraints(self, db_session: Session):
        """Test Task model nullable constraints for required fields.

        Required fields (nullable=False or default NOT NULL):
        - task_id (String)
        - start_time (DateTime)
        - status (String)
        - user_id (Integer, FK)
        - workflow_id (Integer, FK)
        """
        from datetime import datetime

        # Test 1: task_id cannot be NULL
        task_no_id = models.Task(
            task_id=None,  # NULL task_id
            start_time=datetime.now(),
            status="running",
            user_id=1,
            workflow_id=1
        )
        db_session.add(task_no_id)
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 2: start_time cannot be NULL
        task_no_start = models.Task(
            task_id="task_001",
            start_time=None,  # NULL start_time
            status="running",
            user_id=1,
            workflow_id=1
        )
        db_session.add(task_no_start)
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 3: status cannot be NULL
        task_no_status = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status=None,  # NULL status
            user_id=1,
            workflow_id=1
        )
        db_session.add(task_no_status)
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 4: user_id cannot be NULL (implicit NOT NULL for FK without nullable=True)
        task_no_user = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=None,  # NULL user_id
            workflow_id=1
        )
        db_session.add(task_no_user)
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

        # Test 5: workflow_id cannot be NULL (implicit NOT NULL for FK without nullable=True)
        task_no_workflow = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=1,
            workflow_id=None  # NULL workflow_id
        )
        db_session.add(task_no_workflow)
        with pytest.raises(IntegrityError):
            db_session.commit()
        db_session.rollback()

    def test_task_multiple_foreign_keys(self, db_session: Session, sample_user: models.User,
                                       sample_workflow: models.Workflow, sample_plugin: models.Plugin):
        """Test Task model with multiple FK relationships (user, workflow, plugin).

        Task has 3 FK relationships:
        - user_id -> users.id (required)
        - workflow_id -> workflows.id (required)
        - plugin_id -> plugins.id (optional, nullable=True)
        """
        from datetime import datetime

        # Arrange & Act - Create task with all 3 FKs
        task = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id,
            plugin_id=sample_plugin.id  # Optional FK
        )
        db_session.add(task)
        db_session.commit()
        db_session.refresh(task)

        # Assert - All 3 relationships loaded correctly
        assert task.user == sample_user, "task.user should reference sample_user"
        assert task.workflows == sample_workflow, "task.workflows should reference sample_workflow"
        assert task.plugin == sample_plugin, "task.plugin should reference sample_plugin"

        # Assert - Bidirectional relationships
        assert task in sample_user.tasks, "task should be in sample_user.tasks"
        assert task in sample_workflow.tasks, "task should be in sample_workflow.tasks"
        assert task in sample_plugin.tasks, "task should be in sample_plugin.tasks"

        # Test invalid user_id FK constraint
        invalid_task_user = models.Task(
            task_id="task_002",
            start_time=datetime.now(),
            status="running",
            user_id=99999,  # Invalid FK
            workflow_id=sample_workflow.id
        )
        db_session.add(invalid_task_user)
        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()
        assert "foreign key" in str(exc_info.value).lower()
        db_session.rollback()

        # Test invalid workflow_id FK constraint
        invalid_task_workflow = models.Task(
            task_id="task_003",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=99999  # Invalid FK
        )
        db_session.add(invalid_task_workflow)
        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()
        assert "foreign key" in str(exc_info.value).lower()
        db_session.rollback()

        # Test invalid plugin_id FK constraint (optional FK but still enforced)
        invalid_task_plugin = models.Task(
            task_id="task_004",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id,
            plugin_id=99999  # Invalid FK
        )
        db_session.add(invalid_task_plugin)
        with pytest.raises(IntegrityError) as exc_info:
            db_session.commit()
        assert "foreign key" in str(exc_info.value).lower()
        db_session.rollback()

    def test_task_optional_fields(self, db_session: Session, sample_user: models.User,
                                  sample_workflow: models.Workflow):
        """Test Task model optional fields can be NULL.

        Optional fields (nullable=True):
        - end_time (DateTime)
        - algorithm_id (String)
        - plugin_name (String)
        - plugin_id (Integer, FK)
        - task_type (String)
        - plugin_image_uri (String)
        """
        from datetime import datetime

        # Arrange & Act - Create task with only required fields (all optional fields NULL)
        task_minimal = models.Task(
            task_id="task_minimal",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id
            # All optional fields omitted (NULL)
        )
        db_session.add(task_minimal)
        db_session.commit()
        db_session.refresh(task_minimal)

        # Assert - Optional fields are NULL
        assert task_minimal.end_time is None, "end_time should be NULL"
        assert task_minimal.algorithm_id is None, "algorithm_id should be NULL"
        assert task_minimal.plugin_name is None, "plugin_name should be NULL"
        assert task_minimal.plugin_id is None, "plugin_id should be NULL"
        assert task_minimal.task_type is None, "task_type should be NULL"
        assert task_minimal.plugin_image_uri is None, "plugin_image_uri should be NULL"

        # Arrange & Act - Create task with all optional fields set
        task_complete = models.Task(
            task_id="task_complete",
            start_time=datetime.now(),
            end_time=datetime.now(),  # Optional
            status="completed",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id,
            algorithm_id="TENET_v1.0",  # Optional
            plugin_name="TENET",  # Optional (backward compatibility)
            task_type="compile",  # Optional
            plugin_image_uri="docker.io/cellcraft/tenet:1.0"  # Optional
        )
        db_session.add(task_complete)
        db_session.commit()
        db_session.refresh(task_complete)

        # Assert - Optional fields are set correctly
        assert task_complete.end_time is not None, "end_time should be set"
        assert task_complete.algorithm_id == "TENET_v1.0"
        assert task_complete.plugin_name == "TENET"
        assert task_complete.task_type == "compile"
        assert task_complete.plugin_image_uri == "docker.io/cellcraft/tenet:1.0"


# ==============================================================================
# Cascade Behavior Tests
# ==============================================================================

@pytest.mark.unit
@pytest.mark.models
class TestCascadeBehavior:
    """Test cascade deletion and relationship cleanup behavior.

    IMPORTANT: These tests document ACTUAL database behavior, not expected behavior.
    Tests verify current database FK constraint behavior for regression detection.

    Database behavior: ON DELETE SET NULL for FK constraints
    - Parent deletion succeeds without FK constraint violation
    - Child records remain in database with NULL FK values (not deleted)
    - FK columns become nullable despite SQLAlchemy model definitions not having nullable=True

    This behavior creates orphaned records with NULL foreign keys.
    """

    def test_cascade_user_deletion_workflows(self, db_session: Session, sample_user: models.User):
        """Document actual behavior when User with Workflows is deleted.

        Current behavior: User deletion succeeds, workflow.user_id becomes NULL.
        FK constraint: workflow.user_id -> users.id (ON DELETE SET NULL)
        Database behavior: User is deleted, workflows remain with NULL user_id

        Note: workflow.user_id becomes nullable at database level despite model not having nullable=True.
        This is due to Alembic migration or manual database schema configuration.
        """
        # Arrange - Create workflow for user
        workflow = models.Workflow(
            title="Test Workflow",
            workflow_info={"drawflow": {"Home": {"data": {}}}},
            user_id=sample_user.id
        )
        db_session.add(workflow)
        db_session.commit()
        workflow_id = workflow.id
        user_id = sample_user.id

        # Verify workflow exists before deletion
        assert db_session.query(models.Workflow).filter_by(id=workflow_id).first() is not None

        # Act - Delete user
        db_session.delete(sample_user)
        db_session.commit()

        # Clear session to force fresh DB query
        db_session.expire_all()

        # Assert - Document actual behavior
        # User is deleted
        user_after = db_session.query(models.User).filter_by(id=user_id).first()
        assert user_after is None, "User should be deleted"

        # Workflow still exists with NULL user_id (SET NULL behavior)
        workflow_after = db_session.query(models.Workflow).filter_by(id=workflow_id).first()
        assert workflow_after is not None, "ACTUAL BEHAVIOR: Workflow remains after User deletion"
        assert workflow_after.user_id is None, "ACTUAL BEHAVIOR: user_id is SET TO NULL (ON DELETE SET NULL)"

    def test_cascade_user_deletion_files(self, db_session: Session, sample_user: models.User):
        """Document actual behavior when User with Files is deleted.

        Current behavior: User deletion succeeds, file.user_id becomes NULL.
        FK constraint: file.user_id -> users.id (ON DELETE SET NULL)
        Database behavior: User is deleted, files remain with NULL user_id

        Note: file.user_id becomes nullable at database level despite model not having nullable=True.
        """
        # Arrange - Create file for user
        file_record = models.File(
            file_name="test.h5ad",
            file_size="1024",
            file_path="/path/to/test.h5ad",
            folder="/path/to",
            user_id=sample_user.id
        )
        db_session.add(file_record)
        db_session.commit()
        file_id = file_record.id
        user_id = sample_user.id

        # Verify file exists before deletion
        assert db_session.query(models.File).filter_by(id=file_id).first() is not None

        # Act - Delete user
        db_session.delete(sample_user)
        db_session.commit()

        # Clear session to force fresh DB query
        db_session.expire_all()

        # Assert - Document actual behavior
        user_after = db_session.query(models.User).filter_by(id=user_id).first()
        assert user_after is None, "User should be deleted"

        file_after = db_session.query(models.File).filter_by(id=file_id).first()
        assert file_after is not None, "ACTUAL BEHAVIOR: File remains after User deletion"
        assert file_after.user_id is None, "ACTUAL BEHAVIOR: user_id is SET TO NULL (ON DELETE SET NULL)"

    def test_cascade_workflow_deletion_tasks(self, db_session: Session, sample_user: models.User,
                                            sample_workflow: models.Workflow):
        """Document actual behavior when Workflow with Tasks is deleted.

        Current behavior: Workflow deletion succeeds, task.workflow_id is SET TO NULL.
        FK constraint: task.workflow_id -> workflows.id (ON DELETE SET NULL)
        Database behavior: Workflow is deleted, tasks remain with NULL workflow_id

        Note: This is expected SET NULL behavior for optional FK relationships.
        """
        from datetime import datetime

        # Arrange - Create task for workflow
        task = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id
        )
        db_session.add(task)
        db_session.commit()
        task_id = task.id
        workflow_id = sample_workflow.id

        # Verify task exists before deletion
        assert db_session.query(models.Task).filter_by(id=task_id).first() is not None

        # Act - Delete workflow
        db_session.delete(sample_workflow)
        db_session.commit()

        # Clear session to force fresh DB query
        db_session.expire_all()

        # Assert - Document actual behavior
        workflow_after = db_session.query(models.Workflow).filter_by(id=workflow_id).first()
        assert workflow_after is None, "Workflow should be deleted"

        task_after = db_session.query(models.Task).filter_by(id=task_id).first()
        assert task_after is not None, "ACTUAL BEHAVIOR: Task remains after Workflow deletion"
        assert task_after.workflow_id is None, "ACTUAL BEHAVIOR: workflow_id is SET TO NULL (ON DELETE SET NULL)"

    def test_cascade_plugin_deletion_tasks(self, db_session: Session, sample_user: models.User,
                                          sample_workflow: models.Workflow, sample_plugin: models.Plugin):
        """Document actual behavior when Plugin with Tasks is deleted.

        Current behavior: Plugin deletion succeeds, task.plugin_id is SET TO NULL.
        FK constraint: task.plugin_id -> plugins.id (ON DELETE SET NULL, nullable FK)
        Database behavior: Plugin is deleted, tasks remain with NULL plugin_id

        Note: This is expected SET NULL behavior for optional FK relationships.
        plugin_id is already nullable, so SET NULL is appropriate.
        """
        from datetime import datetime

        # Arrange - Create task referencing plugin
        task = models.Task(
            task_id="task_001",
            start_time=datetime.now(),
            status="running",
            user_id=sample_user.id,
            workflow_id=sample_workflow.id,
            plugin_id=sample_plugin.id
        )
        db_session.add(task)
        db_session.commit()
        task_id = task.id
        plugin_id = sample_plugin.id

        # Verify task exists before deletion
        assert db_session.query(models.Task).filter_by(id=task_id).first() is not None

        # Act - Delete plugin
        db_session.delete(sample_plugin)
        db_session.commit()

        # Clear session to force fresh DB query
        db_session.expire_all()

        # Assert - Document actual behavior
        plugin_after = db_session.query(models.Plugin).filter_by(id=plugin_id).first()
        assert plugin_after is None, "Plugin should be deleted"

        task_after = db_session.query(models.Task).filter_by(id=task_id).first()
        assert task_after is not None, "ACTUAL BEHAVIOR: Task remains after Plugin deletion"
        assert task_after.plugin_id is None, "ACTUAL BEHAVIOR: plugin_id is SET TO NULL (ON DELETE SET NULL)"

    def test_many_to_many_association_cleanup(self, db_session: Session, user_factory):
        """Test user_plugin_association table cleanup on User or Plugin deletion.

        Current behavior: Association records CASCADE DELETE (automatically removed).
        Association table: user_plugin_association (user_id, plugin_id)
        FK constraints: Both FKs use CASCADE DELETE behavior

        Database behavior: Deleting User or Plugin automatically removes association records.
        """
        # Arrange - Create user and plugin with association
        user = user_factory(username="test_user", email="test@example.com")

        plugin = models.Plugin(
            name="TestPlugin",
            description="Test plugin for association",
            author="test_author",
            plugin_path="./plugin/local/TestPlugin/",
            plugin_type=None,
            dependencies={"requirements.txt": "numpy==1.24.0"},
            drawflow={"drawflow": {"Home": {"data": {}}}},
            rules={"Snakefile": "rule test:\\n    output: 'test.txt'"},
            use_gpu=False,
            source="local",
            is_editable=True
        )
        db_session.add(plugin)
        db_session.commit()
        db_session.refresh(plugin)

        # Create association through relationship
        user.plugins.append(plugin)
        db_session.commit()

        # Verify association exists
        assert plugin in user.plugins, "plugin should be in user.plugins"
        assert user in plugin.users, "user should be in plugin.users"

        user_id = user.id
        plugin_id = plugin.id

        # Test 1: Delete User and verify association is automatically removed (CASCADE)
        db_session.delete(user)
        db_session.commit()

        # Query association table directly
        from sqlalchemy import text
        result = db_session.execute(
            text("SELECT COUNT(*) FROM user_plugin_association WHERE user_id = :user_id"),
            {"user_id": user_id}
        ).scalar()
        assert result == 0, "Association records should be deleted when User is deleted (CASCADE)"

        # Test 2: Create new user association and test plugin deletion cascade
        user2 = user_factory(username="test_user2", email="test2@example.com")
        user2.plugins.append(plugin)
        db_session.commit()

        # Verify new association exists
        result = db_session.execute(
            text("SELECT COUNT(*) FROM user_plugin_association WHERE plugin_id = :plugin_id"),
            {"plugin_id": plugin_id}
        ).scalar()
        assert result == 1, "New association should exist"

        # Delete plugin and verify association is automatically removed (CASCADE)
        db_session.delete(plugin)
        db_session.commit()

        result = db_session.execute(
            text("SELECT COUNT(*) FROM user_plugin_association WHERE plugin_id = :plugin_id"),
            {"plugin_id": plugin_id}
        ).scalar()
        assert result == 0, "Association records should be deleted when Plugin is deleted (CASCADE)"
