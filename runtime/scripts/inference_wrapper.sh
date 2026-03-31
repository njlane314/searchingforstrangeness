#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RUNTIME_DIR="$(cd "$THIS_DIR/.." && pwd -P)"
INFER_DIR="$RUNTIME_DIR/inference"

if [[ ! -f "$INFER_DIR/infer_bin.py" && -f "$RUNTIME_DIR/infer_bin.py" ]]; then
  INFER_DIR="$RUNTIME_DIR"
fi

if [[ ! -f "$INFER_DIR/infer_bin.py" ]]; then
  echo "ERROR: infer_bin.py not found in RUNTIME_DIR=$RUNTIME_DIR" >&2
  exit 1
fi

cd "$INFER_DIR"
exec python3 -u infer_bin.py "$@"
