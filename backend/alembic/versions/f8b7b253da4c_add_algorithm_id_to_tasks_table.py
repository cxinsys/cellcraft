"""Add algorithm_id to tasks table

Revision ID: f8b7b253da4c
Revises: 
Create Date: 2025-06-16 17:10:46.760881

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = 'f8b7b253da4c'
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('celery_tasksetmeta')
    op.drop_table('celery_taskmeta')
    op.add_column('tasks', sa.Column('algorithm_id', sa.String(), nullable=True))
    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('tasks', 'algorithm_id')
    op.create_table('celery_taskmeta',
    sa.Column('id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('task_id', sa.VARCHAR(length=155), autoincrement=False, nullable=True),
    sa.Column('status', sa.VARCHAR(length=50), autoincrement=False, nullable=True),
    sa.Column('result', postgresql.BYTEA(), autoincrement=False, nullable=True),
    sa.Column('date_done', postgresql.TIMESTAMP(), autoincrement=False, nullable=True),
    sa.Column('traceback', sa.TEXT(), autoincrement=False, nullable=True),
    sa.Column('name', sa.VARCHAR(length=155), autoincrement=False, nullable=True),
    sa.Column('args', postgresql.BYTEA(), autoincrement=False, nullable=True),
    sa.Column('kwargs', postgresql.BYTEA(), autoincrement=False, nullable=True),
    sa.Column('worker', sa.VARCHAR(length=155), autoincrement=False, nullable=True),
    sa.Column('retries', sa.INTEGER(), autoincrement=False, nullable=True),
    sa.Column('queue', sa.VARCHAR(length=155), autoincrement=False, nullable=True),
    sa.PrimaryKeyConstraint('id', name='celery_taskmeta_pkey'),
    sa.UniqueConstraint('task_id', name='celery_taskmeta_task_id_key')
    )
    op.create_table('celery_tasksetmeta',
    sa.Column('id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('taskset_id', sa.VARCHAR(length=155), autoincrement=False, nullable=True),
    sa.Column('result', postgresql.BYTEA(), autoincrement=False, nullable=True),
    sa.Column('date_done', postgresql.TIMESTAMP(), autoincrement=False, nullable=True),
    sa.PrimaryKeyConstraint('id', name='celery_tasksetmeta_pkey'),
    sa.UniqueConstraint('taskset_id', name='celery_tasksetmeta_taskset_id_key')
    )
    # ### end Alembic commands ### 