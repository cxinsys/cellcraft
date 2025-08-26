"""add plugin name source unique constraint

Revision ID: a8f3c2b1d9e4
Revises: ff8810e9f544
Create Date: 2025-08-14 10:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy import text

# revision identifiers, used by Alembic.
revision = 'a8f3c2b1d9e4'
down_revision = 'ff8810e9f544'
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Handle existing data conflicts before adding constraint
    # First, check for duplicates and rename local plugins that conflict with official ones
    connection = op.get_bind()
    
    # Find all plugins with duplicate names but different sources
    result = connection.execute(text("""
        SELECT p1.id, p1.name, p1.source
        FROM plugins p1
        WHERE EXISTS (
            SELECT 1 FROM plugins p2
            WHERE p1.name = p2.name
            AND p1.source != p2.source
            AND p1.id != p2.id
        )
        ORDER BY p1.name, p1.source
    """))
    
    conflicts = result.fetchall()
    
    # Process conflicts - append "_local" to local plugins that conflict with official ones
    processed_names = set()
    for plugin_id, name, source in conflicts:
        if source == 'local' and name not in processed_names:
            # Check if there's an official plugin with the same name
            official_exists = connection.execute(text("""
                SELECT COUNT(*) FROM plugins 
                WHERE name = :name AND source = 'official'
            """), {"name": name}).scalar() > 0
            
            if official_exists:
                # Rename the local plugin
                new_name = f"{name}_local"
                connection.execute(text("""
                    UPDATE plugins 
                    SET name = :new_name 
                    WHERE id = :plugin_id
                """), {"new_name": new_name, "plugin_id": plugin_id})
                processed_names.add(name)
    
    # Check if there's an existing unique constraint on just the name column
    # and drop it if it exists
    inspector = sa.inspect(connection)
    existing_constraints = inspector.get_unique_constraints('plugins')
    
    for constraint in existing_constraints:
        if 'name' in constraint['column_names'] and len(constraint['column_names']) == 1:
            op.drop_constraint(constraint['name'], 'plugins', type_='unique')
    
    # Add the composite unique constraint on (name, source)
    op.create_unique_constraint(
        'uq_plugin_name_source',
        'plugins',
        ['name', 'source']
    )


def downgrade() -> None:
    # Drop the composite unique constraint
    op.drop_constraint('uq_plugin_name_source', 'plugins', type_='unique')
    
    # Note: We don't restore the original names of renamed plugins
    # as that could cause conflicts again
    
    # Optionally, you could add back a unique constraint on just name
    # but this might fail if there are duplicates
    # op.create_unique_constraint('uq_plugin_name', 'plugins', ['name'])