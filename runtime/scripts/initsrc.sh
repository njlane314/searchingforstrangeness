#!/usr/bin/env bash
: "${RUNTIME_BASE_DIR:=}"
: "${ASSETS_BASE_DIR:=}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
DEFAULT_RUNTIME_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"

if [[ -z "$RUNTIME_BASE_DIR" && -n "$ASSETS_BASE_DIR" ]]; then
  RUNTIME_BASE_DIR="$ASSETS_BASE_DIR"
fi

if [[ -z "$RUNTIME_BASE_DIR" && -n "${CONDOR_DIR_INPUT:-}" && -d "${CONDOR_DIR_INPUT}/strangeness/runtime" ]]; then
  RUNTIME_BASE_DIR="${CONDOR_DIR_INPUT}/strangeness/runtime"
elif [[ -z "$RUNTIME_BASE_DIR" && -n "${CONDOR_DIR_INPUT:-}" && -d "${CONDOR_DIR_INPUT}/strangeness/assets" ]]; then
  RUNTIME_BASE_DIR="${CONDOR_DIR_INPUT}/strangeness/assets"
elif [[ -z "$RUNTIME_BASE_DIR" && -d "$PWD/runtime" ]]; then
  RUNTIME_BASE_DIR="$(cd "$PWD/runtime" && pwd -P)"
elif [[ -z "$RUNTIME_BASE_DIR" && -d "$PWD/assets" ]]; then
  RUNTIME_BASE_DIR="$(cd "$PWD/assets" && pwd -P)"
elif [[ -z "$RUNTIME_BASE_DIR" ]]; then
  RUNTIME_BASE_DIR="$DEFAULT_RUNTIME_DIR"
fi

if [[ ! -d "$RUNTIME_BASE_DIR" ]]; then
  echo "ERROR: Could not locate runtime dir. Set RUNTIME_BASE_DIR or ASSETS_BASE_DIR." >&2
  return 1 2>/dev/null || exit 1
fi

if [[ "$PWD" != "$RUNTIME_BASE_DIR" && ! -e "$PWD/runtime" ]]; then
  ln -s "$RUNTIME_BASE_DIR" "$PWD/runtime"
fi

if [[ "$PWD" != "$RUNTIME_BASE_DIR" && ! -e "$PWD/assets" ]]; then
  ln -s "$RUNTIME_BASE_DIR" "$PWD/assets"
fi

export RUNTIME_BASE_DIR
export ASSETS_BASE_DIR="${ASSETS_BASE_DIR:-$RUNTIME_BASE_DIR}"
export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$RUNTIME_BASE_DIR/models}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$RUNTIME_BASE_DIR/detector/badchannels.txt}"
export IA_INFERENCE_WRAPPER="${IA_INFERENCE_WRAPPER:-$RUNTIME_BASE_DIR/scripts/inference_wrapper.sh}"

echo
echo RUNTIME_BASE_DIR=$RUNTIME_BASE_DIR
echo ASSETS_BASE_DIR=$ASSETS_BASE_DIR
echo WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR
echo IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER
echo IA_BADCHANNELS=$IA_BADCHANNELS
echo
