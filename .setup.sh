#!/bin/bash

if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

if ! command -v mrb >/dev/null 2>&1; then
  setup mrb
fi

unsetup_all || true

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_82 -q e17:prof

source /cvmfs/uboone.opensciencegrid.org/products/setup
setup sam_web_client

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbslp

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

UBSIM_EVENTWEIGHT_BASE="${WRK_DIR}/../ubsim/ubsim/EventWeight"

prepend_fhicl_path "${WRK_DIR}/dev"
prepend_fhicl_path "${WRK_DIR}/job"

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
