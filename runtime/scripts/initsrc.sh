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
UBSIM_BASE="${REPO_BASE_DIR}/../ubsim/ubsim"
UBSIM_EVENTWEIGHT_BASE="${REPO_BASE_DIR}/../ubsim/ubsim/EventWeight"
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

  if [[ -n "${FHICL_FILE_PATH:-}" ]]; then
    export FHICL_FILE_PATH="${dir}:${FHICL_FILE_PATH}"
  else
    export FHICL_FILE_PATH="${dir}"
  fi
}

prepend_env_path() {
  local var_name="$1"
  local dir="$2"
  local current_value

  if [[ ! -d "$dir" ]]; then
    return
  fi

  current_value="${!var_name:-}"
  case ":${current_value}:" in
    *:"${dir}":*)
      return
      ;;
  esac

  if [[ -n "${current_value}" ]]; then
    export "${var_name}=${dir}:${current_value}"
  else
    export "${var_name}=${dir}"
  fi
}

prepend_local_product_runtime() {
  local product="$1"
  local product_root="${MRB_INSTALL:-}/${product}"
  local job_dir
  local lib_dir

  if [[ ! -d "${product_root}" ]]; then
    return
  fi

  for job_dir in "${product_root}"/*/job; do
    prepend_fhicl_path "${job_dir}"
  done

  for lib_dir in "${product_root}"/*/*/lib; do
    prepend_env_path CET_PLUGIN_PATH "${lib_dir}"
    prepend_env_path LD_LIBRARY_PATH "${lib_dir}"
  done
}

prepend_fhicl_path "${REPO_BASE_DIR}/dev"
prepend_fhicl_path "${REPO_BASE_DIR}/dev/flux"
prepend_fhicl_path "${REPO_BASE_DIR}/job"
prepend_fhicl_path "${UBSIM_BASE}/Simulation"
prepend_local_product_runtime "ubsim"
prepend_local_product_runtime "ubana"
prepend_fhicl_path "${REDK2NU_BASE}"

if [[ -d "${UBSIM_EVENTWEIGHT_BASE}" ]]; then
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/App"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/numi"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/numi/fluxreader"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/reint"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/genie"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/splines"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/xs"
fi

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
