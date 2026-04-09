#!/usr/bin/env bash

: "${RUNTIME_BASE_DIR:=}"
INITSRC_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
DEFAULT_RUNTIME_DIR="$(cd "$INITSRC_SCRIPT_DIR/.." && pwd -P)"

if [[ -z "$RUNTIME_BASE_DIR" && -n "${CONDOR_DIR_INPUT:-}" && -d "${CONDOR_DIR_INPUT}/strangeness/runtime" ]]; then
  RUNTIME_BASE_DIR="${CONDOR_DIR_INPUT}/strangeness/runtime"
elif [[ -z "$RUNTIME_BASE_DIR" && -d "$PWD/runtime" ]]; then
  RUNTIME_BASE_DIR="$(cd "$PWD/runtime" && pwd -P)"
elif [[ -z "$RUNTIME_BASE_DIR" ]]; then
  RUNTIME_BASE_DIR="$DEFAULT_RUNTIME_DIR"
fi

if [[ ! -d "$RUNTIME_BASE_DIR" ]]; then
  echo "ERROR: Could not locate runtime dir. Set RUNTIME_BASE_DIR." >&2
  return 1 2>/dev/null || exit 1
fi

REPO_BASE_DIR="$(cd "$INITSRC_SCRIPT_DIR/../.." && pwd -P)"
REDK2NU_BASE="${REPO_BASE_DIR}/../redk2nu"

prepend_fhicl_path() {
  local dir="$1"

  if [[ ! -d "$dir" ]]; then
    return
  fi

  case ":${FHICL_FILE_PATH:-}:" in
    *:"${dir}":*)
      return
      ;;
  esac

  export FHICL_FILE_PATH="${dir}${FHICL_FILE_PATH:+:${FHICL_FILE_PATH}}"
}

cleaned_fhicl_path=""
IFS=':' read -r -a FHICL_DIRS <<< "${FHICL_FILE_PATH:-}"
for dir in "${FHICL_DIRS[@]}"; do
  if [[ -z "$dir" ]]; then
    continue
  fi

  case "$dir" in
    "${REPO_BASE_DIR}/job"|\
    "${REPO_BASE_DIR}/job/flux"|\
    "${REPO_BASE_DIR}/job/reinteractions"|\
    */build*/ubana/job|\
    */localProducts*/ubana/*/job)
      continue
      ;;
  esac

  cleaned_fhicl_path="${cleaned_fhicl_path:+${cleaned_fhicl_path}:}${dir}"
done
export FHICL_FILE_PATH="$cleaned_fhicl_path"

prepend_fhicl_path "${REPO_BASE_DIR}/job"
prepend_fhicl_path "${REPO_BASE_DIR}/job/dev"
prepend_fhicl_path "${REPO_BASE_DIR}/job/dev/flux"
prepend_fhicl_path "${REDK2NU_BASE}"

if [[ "$PWD" != "$RUNTIME_BASE_DIR" && ! -e "$PWD/runtime" ]]; then
  ln -s "$RUNTIME_BASE_DIR" "$PWD/runtime"
fi

export RUNTIME_BASE_DIR
export WEIGHTS_BASE_DIR="${WEIGHTS_BASE_DIR:-$RUNTIME_BASE_DIR/models}"
export IA_BADCHANNELS="${IA_BADCHANNELS:-$RUNTIME_BASE_DIR/detector/badchannels.txt}"
export IA_INFERENCE_WRAPPER="${IA_INFERENCE_WRAPPER:-$RUNTIME_BASE_DIR/scripts/inference_wrapper.sh}"

echo
echo RUNTIME_BASE_DIR=$RUNTIME_BASE_DIR
echo WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR
echo IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER
echo IA_BADCHANNELS=$IA_BADCHANNELS
echo FHICL_FILE_PATH=${FHICL_FILE_PATH:-}
echo
