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
import posixpath
import re
from typing import Any, Dict, Iterable, Optional, Set

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


def validate_relative_output_path(value: Any, field: str = "output file name") -> str:
    """결과 디렉터리 기준 상대 경로(하위 디렉터리 허용, ``..``·절대 경로 금지)를 정규화해 반환한다.

    로컬 플러그인 출력명은 작성자가 자유롭게 정할 수 있어 ``plots/out.png`` 같은 하위 경로가 가능하다.
    """
    if not isinstance(value, str) or not value or len(value) > MAX_FILENAME_LENGTH:
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    if value.startswith("/") or re.search(r'[\\"\x00-\x1f\x7f]', value):
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    normalized = posixpath.normpath(value)
    if normalized in (".", "..") or normalized.startswith("../") or "/../" in f"/{normalized}/":
        raise ValidationFailedError(f"Invalid {field}: {value!r}")
    return normalized


def sanitize_filename_part(value: Any) -> str:
    """결과 파일명에 들어가는 파라미터 값에서 경로·따옴표·제어 문자를 ``_``로 바꾼다."""
    return _FORBIDDEN_FILENAME_CHARS.sub("_", str(value))


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


FILE_PARAMETER_TYPES = ("inputFile", "optionalInputFile", "outputFile")
INPUT_FILE_PARAMETER_TYPES = ("inputFile", "optionalInputFile")


def load_plugin_placeholders(plugin_path: str, rule_name: Optional[str] = None) -> Optional[Dict[str, Set[str]]]:
    """플러그인 metadata.json에 선언된 Snakefile placeholder 목록 (읽을 수 없거나 형식이 다르면 None).

    Returns:
        ``{"params": 파일이 아닌 파라미터 이름, "input_files": 입력 파일 파라미터의 defaultValue}``
    """
    try:
        with open(os.path.join(plugin_path, "metadata.json"), "r", encoding="utf-8") as f:
            rules = json.load(f).get("rules", {})
    except (OSError, ValueError, AttributeError):
        return None
    if isinstance(rules, dict):
        rules = list(rules.values())
    if not isinstance(rules, list):
        return None

    placeholders = {"params": set(), "input_files": set()}
    for rule in rules:
        if not isinstance(rule, dict) or (rule_name is not None and rule.get("name") != rule_name):
            continue
        parameters = rule.get("parameters", [])
        if not isinstance(parameters, list):
            continue
        for parameter in parameters:
            if not isinstance(parameter, dict):
                continue
            if parameter.get("type") in INPUT_FILE_PARAMETER_TYPES:
                if isinstance(parameter.get("defaultValue"), str):
                    placeholders["input_files"].add(parameter["defaultValue"])
            elif parameter.get("type") not in FILE_PARAMETER_TYPES and parameter.get("name") is not None:
                placeholders["params"].add(parameter["name"])
    return placeholders


def placeholder_names(placeholders: Optional[Dict[str, Set[str]]], kind: str) -> Optional[Set[str]]:
    return None if placeholders is None else placeholders[kind]


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
