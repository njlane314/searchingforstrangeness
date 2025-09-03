#!/usr/bin/env bash
set -euo pipefail

ASSETS_BASE_DIR="${ASSETS_BASE_DIR:-}"
if [[ -z "${ASSETS_BASE_DIR}" ]]; then
  if [[ -d "$PWD/assets" ]]; then
    ASSETS_BASE_DIR="$(cd "$PWD/assets" && pwd)"
  else
    THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CAND="$(cd "$THIS_DIR/.." && pwd)"
    if [[ -d "$CAND/calib" && -d "$CAND/weights" ]]; then
      ASSETS_BASE_DIR="$CAND"
    else
      echo "ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR." >&2
      exit 1
    fi
  fi
fi
export ASSETS_BASE_DIR

export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$ASSETS_BASE_DIR/weights}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$ASSETS_BASE_DIR/calib/badchannels.txt}"
export IA_INFERENCE_WRAPPER="${IA_INFERENCE_WRAPPER:-$ASSETS_BASE_DIR/scripts/run_strangeness_inference.sh}"

export PYTHONPATH="${PYTHONPATH:-}:$ASSETS_BASE_DIR:$ASSETS_BASE_DIR/models:$ASSETS_BASE_DIR/scripts"

if [[ -n "${APPTAINER_BINDPATH:-}" ]]; then
  export APPTAINER_BINDPATH="$ASSETS_BASE_DIR,${APPTAINER_BINDPATH}"
else
  export APPTAINER_BINDPATH="$ASSETS_BASE_DIR"
fi

if [[ -n "${FW_SEARCH_PATH:-}" ]]; then
  export FW_SEARCH_PATH="$ASSETS_BASE_DIR:$FW_SEARCH_PATH"
else
  export FW_SEARCH_PATH="$ASSETS_BASE_DIR"
fi

echo "[init] ASSETS_BASE_DIR=$ASSETS_BASE_DIR"
echo "[init] WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR"
echo "[init] IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER"
echo "[init] IA_BADCHANNELS=$IA_BADCHANNELS"
echo "[init] PYTHONPATH=$PYTHONPATH"
echo "[init] APPTAINER_BINDPATH=$APPTAINER_BINDPATH"

