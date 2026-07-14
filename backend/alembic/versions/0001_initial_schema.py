"""initial schema - create all tables

Revision ID: 0001_initial
Revises:
Create Date: 2026-03-24 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = '0001_initial'
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # === 1. Create enum types ===
    plugintype_enum = sa.Enum('ANALYSIS', 'VISUALIZATION', name='plugintype')
    plugintype_enum.create(op.get_bind(), checkfirst=True)

    # === 2. Create tables (order: independent tables first, then dependent tables) ===

    # --- users (no FK dependencies) ---
    op.create_table(
        'users',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('username', sa.String(), nullable=False),
        sa.Column('email', sa.String(), nullable=False),
        sa.Column('hashed_password', sa.String(), nullable=False),
        sa.Column('is_active', sa.Boolean(), default=True),
        sa.Column('is_superuser', sa.Boolean(), default=False),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
    )
    op.create_index('ix_users_id', 'users', ['id'])
    op.create_index('ix_users_email', 'users', ['email'], unique=True)

    # --- plugins (no FK dependencies) ---
    op.create_table(
        'plugins',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('description', sa.String(), nullable=False),
        sa.Column('author', sa.String(), nullable=False),
        sa.Column('plugin_path', sa.String(), nullable=False),
        # create_type=False is only honored by postgresql.ENUM (generic sa.Enum
        # ignores it and emits a duplicate CREATE TYPE, breaking fresh installs).
        sa.Column('plugin_type', postgresql.ENUM('ANALYSIS', 'VISUALIZATION', name='plugintype', create_type=False), nullable=True),
        sa.Column('dependencies', postgresql.JSONB(), nullable=True),
        sa.Column('drawflow', postgresql.JSONB(), nullable=False),
        sa.Column('rules', postgresql.JSONB(), nullable=False),
        sa.Column('use_gpu', sa.Boolean(), default=False),
        sa.Column('source', sa.String(), nullable=False, server_default='local'),
        sa.Column('is_editable', sa.Boolean(), nullable=False, server_default='true'),
        sa.Column('version', sa.String(), nullable=True),
        sa.Column('submodule_path', sa.String(), nullable=True),
        sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.UniqueConstraint('name', 'source', name='uq_plugin_name_source'),
    )
    op.create_index('ix_plugins_id', 'plugins', ['id'])

    # --- files (FK: users) ---
    op.create_table(
        'files',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('file_name', sa.String(), nullable=False),
        sa.Column('file_size', sa.String(), nullable=False),
        sa.Column('file_path', sa.String(), nullable=False),
        sa.Column('folder', sa.String(), nullable=False),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
    )
    op.create_index('ix_files_id', 'files', ['id'])

    # --- workflows (FK: users) ---
    op.create_table(
        'workflows',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('title', sa.String(), nullable=False),
        sa.Column('thumbnail', sa.Text(), nullable=True),
        sa.Column('workflow_info', postgresql.JSONB(), nullable=False),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
    )
    op.create_index('ix_workflows_id', 'workflows', ['id'])

    # --- tasks (FK: users, workflows, plugins) ---
    op.create_table(
        'tasks',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('task_id', sa.String(), nullable=False),
        sa.Column('start_time', sa.DateTime(timezone=True), nullable=False),
        sa.Column('end_time', sa.DateTime(timezone=True), nullable=True),
        sa.Column('status', sa.String(), nullable=False),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('workflow_id', sa.Integer(), sa.ForeignKey('workflows.id'), nullable=True),
        sa.Column('algorithm_id', sa.String(), nullable=True),
        sa.Column('plugin_name', sa.String(), nullable=True),
        sa.Column('plugin_id', sa.Integer(), sa.ForeignKey('plugins.id'), nullable=True),
        sa.Column('task_type', sa.String(), nullable=True),
        sa.Column('plugin_image_uri', sa.String(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
    )
    op.create_index('ix_tasks_plugin_id', 'tasks', ['plugin_id'])

    # --- user_plugin_association (FK: users, plugins) ---
    op.create_table(
        'user_plugin_association',
        sa.Column('user_id', sa.Integer(), sa.ForeignKey('users.id'), nullable=True),
        sa.Column('plugin_id', sa.Integer(), sa.ForeignKey('plugins.id'), nullable=True),
    )


def downgrade() -> None:
    # Drop tables in reverse dependency order
    op.drop_table('user_plugin_association')
    op.drop_table('tasks')
    op.drop_table('workflows')
    op.drop_table('files')
    op.drop_table('plugins')
    op.drop_table('users')

    # Drop enum types
    plugintype_enum = sa.Enum(name='plugintype')
    plugintype_enum.drop(op.get_bind(), checkfirst=True)
