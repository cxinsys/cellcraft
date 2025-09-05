"""add plugin_id foreign key to tasks table

Revision ID: add_plugin_id_fk_to_tasks
Revises: a1b2c3d4e5f6
Create Date: 2025-09-02 10:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy import text

# revision identifiers, used by Alembic.
revision = 'add_plugin_id_fk_to_tasks'
down_revision = 'a1b2c3d4e5f6'
branch_labels = None
depends_on = None


def upgrade() -> None:
    """Add plugin_id foreign key to tasks table and migrate existing data"""
    
    # Step 1: Add the plugin_id column as nullable
    op.add_column('tasks', sa.Column('plugin_id', sa.Integer(), nullable=True))
    
    # Step 2: Create the foreign key constraint
    op.create_foreign_key(
        constraint_name='fk_tasks_plugin_id',
        source_table='tasks',
        referent_table='plugins',
        local_cols=['plugin_id'],
        remote_cols=['id'],
        ondelete='SET NULL'  # If plugin is deleted, set plugin_id to NULL
    )
    
    # Step 3: Data migration - populate plugin_id based on plugin_name
    # This query matches tasks to plugins by name, prioritizing exact matches
    connection = op.get_bind()
    connection.execute(text("""
        UPDATE tasks 
        SET plugin_id = plugins.id
        FROM plugins 
        WHERE tasks.plugin_name IS NOT NULL 
        AND tasks.plugin_name = plugins.name
        AND tasks.plugin_id IS NULL;
    """))
    
    # Step 4: Handle case-insensitive matches for any remaining unmatched tasks
    connection.execute(text("""
        UPDATE tasks 
        SET plugin_id = plugins.id
        FROM plugins 
        WHERE tasks.plugin_name IS NOT NULL 
        AND LOWER(tasks.plugin_name) = LOWER(plugins.name)
        AND tasks.plugin_id IS NULL;
    """))
    
    # Step 5: Create index for better query performance
    op.create_index(
        index_name='ix_tasks_plugin_id',
        table_name='tasks',
        columns=['plugin_id']
    )


def downgrade() -> None:
    """Remove plugin_id foreign key from tasks table"""
    
    # Drop the index first
    op.drop_index('ix_tasks_plugin_id', table_name='tasks')
    
    # Drop the foreign key constraint
    op.drop_constraint('fk_tasks_plugin_id', 'tasks', type_='foreignkey')
    
    # Drop the plugin_id column
    op.drop_column('tasks', 'plugin_id')