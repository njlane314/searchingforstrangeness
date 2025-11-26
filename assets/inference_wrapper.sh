#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$THIS_DIR/.." && pwd)"

if [[ -z "${ASSETS_BASE_DIR:-}" ]]; then
  if [[ -d "$ROOT_DIR/assets" ]]; then
    ASSETS_BASE_DIR="$ROOT_DIR/assets"
  else
    ASSETS_BASE_DIR="$ROOT_DIR"
  fi
fi

if [[ ! -d "$ASSETS_BASE_DIR" ]]; then
  echo "ERROR: ASSETS_BASE_DIR invalid: $ASSETS_BASE_DIR" >&2
  exit 1
fi
cd "$ASSETS_BASE_DIR"
exec python3 -u infer_bin.py "$@"
