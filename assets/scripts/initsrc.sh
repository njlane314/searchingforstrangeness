#!/usr/bin/env bash
: "${ASSETS_BASE_DIR:=}"
if [[ -n "${CONDOR_DIR_INPUT:-}" && -d "${CONDOR_DIR_INPUT}/strangeness/assets" ]]; then
  ASSETS_BASE_DIR="${CONDOR_DIR_INPUT}/strangeness/assets"
elif [[ -z "${ASSETS_BASE_DIR}" ]]; then
  if [[ -d "$PWD/assets" ]]; then
    ASSETS_BASE_DIR="$(cd "$PWD/assets" && pwd)"
  else
    THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CAND="$(cd "$THIS_DIR/.." && pwd)"
    if [[ -d "$CAND/calib" && -d "$CAND/weights" ]]; then
      ASSETS_BASE_DIR="$CAND"
    else
      echo "ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR." >&2
      return 1
    fi
  fi
fi
export ASSETS_BASE_DIR
export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$ASSETS_BASE_DIR/weights}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$ASSETS_BASE_DIR/calib/badchannels.txt}"
export IA_INFERENCE_WRAPPER="${IA_INFERENCE_WRAPPER:-$ASSETS_BASE_DIR/scripts/run_strangeness_inference.sh}"
echo
echo ASSETS_BASE_DIR=$ASSETS_BASE_DIR
echo WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR
echo IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER
echo IA_BADCHANNELS=$IA_BADCHANNELS
echo
