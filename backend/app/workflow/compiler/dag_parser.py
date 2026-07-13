"""Backward-compat facade for the Snakefile DAG parser & status tracker (PR-10 split).

The 2,227-line monolith was split into focused modules. The bulk was the
``SnakemakeRuleStatusTracker`` (~1,455 lines), which exceeds the 800-line target on
its own, so it was decomposed by responsibility rather than by the doc's original
parsing/dag/rules/serialize grouping (``SnakemakeDAGParser`` is one cohesive class
already covering parse/dag/serialize):

- ``compiler/cache.py``           — ``DAGCache``/``LogFileCache`` + cache mgmt helpers
- ``compiler/parsing.py``         — ``DAGParsingError`` + ``SnakemakeDAGParser``
- ``compiler/status.py``          — ``SnakemakeRuleStatusTracker`` (core + log IO + post-process)
- ``compiler/status_analysis.py`` — per-rule verdict mixin (composed into the tracker)
- ``compiler/status_progress.py`` — progress/timing/bottleneck mixin (composed into the tracker)

Existing call sites and test patch targets that reference these names on
``app.workflow.compiler.dag_parser`` (e.g. ``dag_parser.SnakemakeDAGParser``) keep
working via the re-exports below. Class/function signatures and behavior are unchanged.
"""
from app.workflow.compiler.cache import (  # noqa: F401
    DAGCache,
    LogFileCache,
    _dag_cache,
    _log_cache,
    clear_all_caches,
    get_cache_stats,
    set_cache_config,
    get_file_hash,
)
from app.workflow.compiler.parsing import (  # noqa: F401
    DAGParsingError,
    SnakemakeDAGParser,
)
from app.workflow.compiler.status import (  # noqa: F401
    SnakemakeRuleStatusTracker,
)
