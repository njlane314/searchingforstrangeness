#!/usr/bin/env bash
set -euo pipefail
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"
SCRIPT="${SELF_DIR}/../pvdv2d_me.py"
exec "$PY" "$SCRIPT" "$@"
