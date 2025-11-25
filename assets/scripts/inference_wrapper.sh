#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSETS_DIR="$(cd "$THIS_DIR/.." && pwd)"

if [[ ! -f "$ASSETS_DIR/infer_bin.py" ]]; then
  echo "ERROR: infer_bin.py not found in ASSETS_DIR=$ASSETS_DIR" >&2
  exit 1
fi

cd "$ASSETS_DIR"
exec python3 -u infer_bin.py "$@"
