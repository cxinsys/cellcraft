"""Central model registry (PR-4, phase-2-domain-skeleton).

Import every domain's models so Alembic's ``env.py`` and ``create_all``-style
callers see the complete metadata on a single ``Base``.

This lives in a dedicated module (not ``db/base.py``) to avoid a circular
import: the domain ``models`` modules import ``Base`` from ``app.db.base``, so
``app.db.base`` cannot import them back at module load time.

Import order matters: ``user.models`` defines ``user_plugin_association``,
which ``plugin.models`` references, so ``user`` must be imported before
``plugin``.
"""
from app.db.base import Base  # noqa: F401 — re-exported for convenience
from app.user.models import User  # noqa: F401
from app.file.models import File  # noqa: F401
from app.workflow.models import Workflow  # noqa: F401
from app.task.models import Task  # noqa: F401
from app.plugin.models import Plugin  # noqa: F401
