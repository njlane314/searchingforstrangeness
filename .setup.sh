#!/bin/bash

if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

unsetup_all || true

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
# Keep the UPS stack aligned with the local ubana/ubsim v08_00_00_82 area.
setup uboonecode v08_00_00_82 -q e17:prof

source /cvmfs/uboone.opensciencegrid.org/products/setup
if [ -z "${MRB_DIR:-}" ] || ! command -v mrb >/dev/null 2>&1; then
  setup mrb
fi
setup sam_web_client

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbsetenv

prepend_fhicl_path() {
  local dir="$1"
  if [ ! -d "${dir}" ]; then
    return
  fi

  case ":${FHICL_FILE_PATH:-}:" in
    *:"${dir}":*)
      return
      ;;
  esac

  if [ -n "${FHICL_FILE_PATH:-}" ]; then
    export FHICL_FILE_PATH="${dir}:${FHICL_FILE_PATH}"
  else
    export FHICL_FILE_PATH="${dir}"
  fi
}

prepend_env_path() {
  local var_name="$1"
  local dir="$2"
  local current_value

  if [ ! -d "${dir}" ]; then
    return
  fi

  current_value="${!var_name:-}"
  case ":${current_value}:" in
    *:"${dir}":*)
      return
      ;;
  esac

  if [ -n "${current_value}" ]; then
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

  if [ ! -d "${product_root}" ]; then
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

UBSIM_BASE="${WRK_DIR}/../ubsim/ubsim"
UBSIM_EVENTWEIGHT_BASE="${WRK_DIR}/../ubsim/ubsim/EventWeight"

prepend_fhicl_path "${WRK_DIR}/dev"
prepend_fhicl_path "${WRK_DIR}/dev/flux"
prepend_fhicl_path "${WRK_DIR}/job"
prepend_fhicl_path "${WRK_DIR}/job/flux"
prepend_fhicl_path "${WRK_DIR}/job/reinteractions"
prepend_fhicl_path "${UBSIM_BASE}/Simulation"
prepend_local_product_runtime "ubsim"
prepend_local_product_runtime "ubana"

if [ -d "${UBSIM_EVENTWEIGHT_BASE}" ]; then
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/App"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/numi"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/numi/fluxreader"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/reint"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/genie"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/splines"
  prepend_fhicl_path "${UBSIM_EVENTWEIGHT_BASE}/jobs/xs"
fi

# For SAM write operations on the GPVM, reset any stale bearer-token env first:
# unset BEARER_TOKEN
# unset BEARER_TOKEN_FILE
# htdestroytoken || true
# htgettoken -a htvaultprod.fnal.gov -i uboone

echo FHICL_FILE_PATH=${FHICL_FILE_PATH}
which jobsub_submit
