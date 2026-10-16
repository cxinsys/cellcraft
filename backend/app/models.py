"""Aggregated model namespace (PR-4, phase-2-domain-skeleton).

Each ORM model now lives in its owning domain (``app.<domain>.models``). This
module re-exports them under a single namespace so existing
``from app import models`` / ``models.User`` call sites keep working without a
per-symbol rewrite. It also pulls in the full registry so ``Base.metadata`` is
complete whenever any caller imports ``models``.
"""
from app.user.models import User, user_plugin_association  # noqa: F401
from app.file.models import File  # noqa: F401
from app.workflow.models import Workflow  # noqa: F401
from app.task.models import Task  # noqa: F401
from app.plugin.models import Plugin  # noqa: F401
