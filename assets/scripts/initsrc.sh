#!/usr/bin/env bash
: "${ASSETS_BASE_DIR:=}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "$ASSETS_BASE_DIR" && -n "${CONDOR_DIR_INPUT:-}" && -d "${CONDOR_DIR_INPUT}/strangeness/assets" ]]; then
  ASSETS_BASE_DIR="${CONDOR_DIR_INPUT}/strangeness/assets"
elif [[ -z "$ASSETS_BASE_DIR" && -d "$PWD/assets" ]]; then
  ASSETS_BASE_DIR="$(cd "$PWD/assets" && pwd)"
elif [[ -z "$ASSETS_BASE_DIR" ]]; then
  ASSETS_BASE_DIR="$ROOT_DIR/assets"
fi

if [[ ! -d "$ASSETS_BASE_DIR" ]]; then
  echo "ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR." >&2
  return 1
fi

if [[ ! -e "$PWD/assets" ]]; then
  ln -s "$ASSETS_BASE_DIR" "$PWD/assets"
fi

export ASSETS_BASE_DIR
export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$ASSETS_BASE_DIR/weights}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$ASSETS_BASE_DIR/badchannels.txt}"
export IA_INFERENCE_WRAPPER="${IA_INFERENCE_WRAPPER:-$ROOT_DIR/app/inference_wrapper.sh}"

echo
echo ASSETS_BASE_DIR=$ASSETS_BASE_DIR
echo WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR
echo IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER
echo IA_BADCHANNELS=$IA_BADCHANNELS
echo
