#!/usr/bin/env bash
# Per-installation credentials for the internal services (PostgreSQL, RabbitMQ, Redis).
#
# Creates or completes <repo>/.env, which docker compose reads automatically:
#   - new install: random passwords are generated
#   - existing install: values already in .env are kept, only missing keys are added
#   - .env missing but the PostgreSQL volume exists: the volume was initialized with
#     the legacy default account, so missing PostgreSQL values fall back to it
#
# Usage:
#   scripts/ensure-env.sh                 create/complete .env
#   source scripts/ensure-env.sh          define the cellcraft_* functions only
#
# Values are never printed. The file is kept at mode 0600.

CELLCRAFT_ROOT="${CELLCRAFT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CELLCRAFT_ENV_FILE="${CELLCRAFT_ENV_FILE:-${CELLCRAFT_ROOT}/.env}"

# 기존 설치(v1.0.7 이하) DB 볼륨이 초기화될 때 사용된 기본 계정
CELLCRAFT_LEGACY_POSTGRES_USER="cellcraft_admin"
CELLCRAFT_LEGACY_POSTGRES_PASSWORD="cellcraft_admin"
CELLCRAFT_LEGACY_POSTGRES_DB="cellcraft"

_cellcraft_log() {
    local level="$1"; shift
    if declare -F "log_${level}" >/dev/null; then
        "log_${level}" "$*"
    else
        echo "[${level}] $*" >&2
    fi
}

# Print the value of KEY from an env file ("KEY=value" or "KEY = value"), empty if unset.
cellcraft_env_get() {
    local key="$1" file="${2:-$CELLCRAFT_ENV_FILE}" value
    [[ -f "$file" ]] || return 0
    value="$(sed -nE "s/^[[:space:]]*(export[[:space:]]+)?${key}[[:space:]]*=[[:space:]]*(.*)$/\2/p" "$file" | tail -n 1)"
    value="${value%$'\r'}"
    value="${value%"${value##*[![:space:]]}"}"
    if [[ "$value" =~ ^\"(.*)\"$ || "$value" =~ ^\'(.*)\'$ ]]; then
        value="${BASH_REMATCH[1]}"
    fi
    printf '%s' "$value"
}

cellcraft_random_secret() {
    od -An -N24 -tx1 /dev/urandom | tr -d ' \n'
}

# docker compose project name (volume prefix): COMPOSE_PROJECT_NAME or the repo directory name
cellcraft_project_name() {
    local name="${COMPOSE_PROJECT_NAME:-$(cellcraft_env_get COMPOSE_PROJECT_NAME)}"
    name="${name:-$(basename "$CELLCRAFT_ROOT")}"
    printf '%s' "$name" | tr '[:upper:]' '[:lower:]' | tr -cd 'a-z0-9_-'
}

cellcraft_volume_exists() {
    docker volume inspect "$(cellcraft_project_name)_$1" >/dev/null 2>&1
}

cellcraft_ensure_env() {
    local file="$CELLCRAFT_ENV_FILE"
    local -a added=()
    local block="" legacy_db=false

    _cellcraft_add() {
        block+="$1=$2"$'\n'
        added+=("$1")
    }

    if [[ -z "$(cellcraft_env_get POSTGRES_PASSWORD)" ]]; then
        if ! docker info >/dev/null 2>&1; then
            _cellcraft_log error "Docker is not running: cannot check for an existing database volume. Start Docker and retry."
            return 1
        fi
        cellcraft_volume_exists postgres_data && legacy_db=true
    fi

    if [[ -z "$(cellcraft_env_get POSTGRES_USER)" ]]; then
        _cellcraft_add POSTGRES_USER "$CELLCRAFT_LEGACY_POSTGRES_USER"
    fi
    if [[ -z "$(cellcraft_env_get POSTGRES_DB)" ]]; then
        _cellcraft_add POSTGRES_DB "$CELLCRAFT_LEGACY_POSTGRES_DB"
    fi
    if [[ -z "$(cellcraft_env_get POSTGRES_PASSWORD)" ]]; then
        if [[ "$legacy_db" == "true" ]]; then
            _cellcraft_add POSTGRES_PASSWORD "$CELLCRAFT_LEGACY_POSTGRES_PASSWORD"
            _cellcraft_log warning "Existing database volume found without .env: using the previous default database account."
        else
            _cellcraft_add POSTGRES_PASSWORD "$(cellcraft_random_secret)"
        fi
    fi
    if [[ -z "$(cellcraft_env_get RABBITMQ_USER)" ]]; then
        _cellcraft_add RABBITMQ_USER "cellcraft"
    fi
    if [[ -z "$(cellcraft_env_get RABBITMQ_PASSWORD)" ]]; then
        _cellcraft_add RABBITMQ_PASSWORD "$(cellcraft_random_secret)"
    fi
    if [[ -z "$(cellcraft_env_get REDIS_PASSWORD)" ]]; then
        _cellcraft_add REDIS_PASSWORD "$(cellcraft_random_secret)"
    fi

    if (( ${#added[@]} > 0 )); then
        (
            umask 077
            if [[ -s "$file" && -n "$(tail -c 1 "$file")" ]]; then
                echo >> "$file"
            fi
            {
                echo "# Added by scripts/ensure-env.sh on $(date '+%Y-%m-%d %H:%M:%S') (per-installation credentials)"
                printf '%s' "$block"
            } >> "$file"
        ) || return 1
        _cellcraft_log info "Updated $(basename "$file"): ${added[*]}"
    fi
    chmod 600 "$file" 2>/dev/null || true
    unset -f _cellcraft_add
}

# Apply RABBITMQ_USER/PASSWORD to a running broker.
# RABBITMQ_DEFAULT_USER only takes effect when the data volume is first created, so an
# upgraded installation still has the old guest account: create the configured user,
# grant it access and remove guest. Safe to run on every start.
cellcraft_sync_rabbitmq_user() {
    local compose_file="$1"
    local user password i
    user="$(cellcraft_env_get RABBITMQ_USER)"
    password="$(cellcraft_env_get RABBITMQ_PASSWORD)"
    [[ -n "$user" && -n "$password" ]] || return 0

    local -a compose=(docker compose -f "$compose_file")
    local -a ctl=("${compose[@]}" exec -T rabbitmq rabbitmqctl -s)

    for i in $(seq 1 60); do
        "${compose[@]}" exec -T rabbitmq rabbitmq-diagnostics -q ping >/dev/null 2>&1 && break
        if (( i == 60 )); then
            _cellcraft_log error "RabbitMQ did not become ready; broker account not synchronized"
            return 1
        fi
        sleep 2
    done

    local users changed=false
    users="$("${ctl[@]}" list_users 2>/dev/null | awk '{print $1}')" || return 1
    if ! grep -qxF -- "$user" <<<"$users"; then
        "${ctl[@]}" add_user "$user" "$password" >/dev/null || return 1
        changed=true
    elif ! "${ctl[@]}" authenticate_user "$user" "$password" >/dev/null 2>&1; then
        "${ctl[@]}" change_password "$user" "$password" >/dev/null || return 1
        changed=true
    fi
    "${ctl[@]}" set_permissions -p / "$user" ".*" ".*" ".*" >/dev/null || return 1
    "${ctl[@]}" set_user_tags "$user" administrator >/dev/null || return 1
    if [[ "$user" != "guest" ]] && grep -qxF guest <<<"$users"; then
        "${ctl[@]}" delete_user guest >/dev/null || return 1
        changed=true
    fi

    if [[ "$changed" == "true" ]]; then
        _cellcraft_log info "RabbitMQ account updated; restarting celery to reconnect"
        "${compose[@]}" restart celery >/dev/null 2>&1 || return 1
    fi
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail
    cellcraft_ensure_env
fi
