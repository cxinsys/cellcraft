"""Server-side validation of client-supplied workflow values (#137).

The frontend sends Drawflow node ids, bare file names and plugin parameters that
the backend turns into filesystem paths (relative to the backend root, e.g.
``./user/<username>/workflow_<id>/algorithm_<node>/results``) and into Snakefile
content. These helpers reject values that could escape the user's directory or
alter the generated Snakefile, while accepting every value the frontend
legitimately produces:

- node ids: Drawflow auto-increment integers (``3`` / ``"3"``)
- file names: bare names without directory parts (uploads are sanitized to
  ``[A-Za-z0-9_.-]``; shared tutorial files keep their on-disk names)
- node data: ``ScatterPlot`` / ``<node id>`` / ``json``
"""
import json
import os
import re
from typing import Any, Iterable, Optional, Set

from app.core.exceptions import NotFoundError, PermissionDeniedError, ValidationFailedError
from app.file.security import validate_file_path

NODE_ID_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
NODE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
NODE_DATA_EXTENSIONS = {"json", "txt", "tsv", "csv"}
MAX_FILENAME_LENGTH = 255
# Snakefile 문자열 리터럴 안에 그대로 들어가므로 따옴표·역슬래시·제어 문자도 허용하지 않는다
_FORBIDDEN_FILENAME_CHARS = re.compile(r'[/\\"\x00-\x1f\x7f]')
# Snakemake 예약 네임스페이스 placeholder ({input.x}, {params.y} ...)
_SNAKEMAKE_NAMESPACES = ("input.", "output.", "params.", "log.", "wildcards.", "resources.", "threads.")


def validate_node_id(value: Any, field: str = "node_id") -> str:
    """Drawflow 노드 ID (정수 또는 숫자 문자열 등 경로 구분자 없는 식별자)."""
    text = str(value) if value is not None else ""
    if isinstance(value, bool) or not NODE_ID_PATTERN.match(text):
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    return text


def validate_node_name(value: Any, field: str = "node_name") -> str:
    text = str(value) if value is not None else ""
    if not NODE_NAME_PATTERN.match(text):
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    return text


def validate_node_data_extension(value: Any) -> str:
    if value not in NODE_DATA_EXTENSIONS:
        raise ValidationFailedError(f"Invalid file_extension: {value!r}")
    return value


def validate_bare_filename(value: Any, field: str = "file name") -> str:
    """디렉터리 구성요소가 없는 파일명인지 확인한다."""
    if not isinstance(value, str) or value in ("", ".", ".."):
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    if len(value) > MAX_FILENAME_LENGTH or _FORBIDDEN_FILENAME_CHARS.search(value):
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    return value


def ensure_within(base_dir: str, path: str) -> str:
    """``path``가 ``base_dir`` 내부인지 확인한다 (symlink 해석 포함)."""
    try:
        validate_file_path(base_dir, path)
    except Exception:
        raise PermissionDeniedError("Access denied: path is outside the allowed directory")
    return path


def resolve_listed_file(directory: str, filename: Any) -> str:
    """디렉터리 목록에 실제로 존재하는 파일명만 경로로 반환한다 (없으면 404)."""
    validate_bare_filename(filename, "filename")
    if not os.path.isdir(directory) or filename not in os.listdir(directory):
        raise NotFoundError(f"File not found: {filename}")
    return os.path.join(directory, filename)


def load_plugin_parameter_names(plugin_path: str, rule_name: Optional[str] = None) -> Optional[Set[str]]:
    """플러그인 metadata.json에 선언된 파라미터 이름 목록 (읽을 수 없으면 None)."""
    try:
        with open(os.path.join(plugin_path, "metadata.json"), "r", encoding="utf-8") as f:
            rules = json.load(f).get("rules", {})
    except (OSError, ValueError, AttributeError):
        return None
    names = set()
    for rule in rules.values():
        if rule_name is not None and rule.get("name") != rule_name:
            continue
        for parameter in rule.get("parameters", []):
            if parameter.get("name") is not None:
                names.add(parameter["name"])
    return names


def filter_plugin_params(params: dict, allowed_names: Optional[Iterable[str]]) -> dict:
    """Snakefile에 치환할 파라미터를 metadata에 선언된 이름으로 제한한다.

    metadata를 읽지 못한 경우에는 Snakemake 예약 네임스페이스 placeholder만 제외한다.
    """
    if allowed_names is not None:
        allowed = set(allowed_names)
        return {name: value for name, value in params.items() if name in allowed}
    return {
        name: value for name, value in params.items()
        if isinstance(name, str) and not name.startswith(_SNAKEMAKE_NAMESPACES)
    }
