#!/usr/bin/env bash
set -euo pipefail

# Resolve ASSETS_BASE_DIR robustly from this script path
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${ASSETS_BASE_DIR:="$(cd "$THIS_DIR/.." && pwd)"}"

# Sanity check the expected subdirs
if [[ ! -d "$ASSETS_BASE_DIR/calib" || ! -d "$ASSETS_BASE_DIR/weights" ]]; then
  echo "ERROR: ASSETS_BASE_DIR invalid: $ASSETS_BASE_DIR" >&2
  exit 1
fi

export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$ASSETS_BASE_DIR/weights}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$ASSETS_BASE_DIR/calib/badchannels.txt}"
if [[ ! -f "$IA_BADCHANNELS" ]]; then
  echo "ERROR: Missing badchannels file at $IA_BADCHANNELS" >&2
  exit 1
fi
export IA_INFERENCE_WRAPPER="$THIS_DIR/run_strangeness_inference.sh"

# Make sure Python will see our packages
export PYTHONPATH="$ASSETS_BASE_DIR:$ASSETS_BASE_DIR/models:$ASSETS_BASE_DIR/scripts${PYTHONPATH:+:$PYTHONPATH}"

# Optional: avoid surprises from a host Python
unset PYTHONHOME

echo "[init] ASSETS_BASE_DIR=$ASSETS_BASE_DIR"
echo "[init] PYTHONPATH=$PYTHONPATH"

cd "$ASSETS_BASE_DIR"
exec python3 -u "$ASSETS_BASE_DIR/scripts/run_inference.py" "$@"
