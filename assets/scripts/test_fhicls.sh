#!/bin/bash

set -u
set -o pipefail

source "assets/scripts/initsrc.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export FHICL_FILE_PATH="${SCRIPT_DIR}/job:${SCRIPT_DIR}:${FHICL_FILE_PATH:-}"

SAM_DEF="prod_strange_resample_fhc_run2_fhc_reco2_reco2"
NUM_EVENTS=1

FHICLS=(
  "${SCRIPT_DIR}/job/imaging_inference_nuselection.fcl"

  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_cv.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_syst_singleknobs.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_extragenieall_1.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_extragenieall_2.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_extragenieall_3.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_extragenieall_4.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_fhc_extragenieall_5.fcl"

  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_cv.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_syst_singleknobs.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_extragenieall_1.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_extragenieall_2.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_extragenieall_3.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_extragenieall_4.fcl"
  "${SCRIPT_DIR}/job/run_eventweight_numi_rhc_extragenieall_5.fcl"
)

FILE="$(samweb list-files defname:"${SAM_DEF}" | head -n 1 || true)"
if [[ -z "${FILE}" ]]; then
  echo "No files returned from SAM definition: ${SAM_DEF}"
  exit 1
fi

FILEDIR="$(samweb locate-file "${FILE}" | grep -o '/pnfs/.*' | sed 's/(\([0-9]*@[a-z0-9]*\))//g' | head -n 1 || true)"
if [[ -z "${FILEDIR}" ]]; then
  echo "Could not locate PNFS directory for file: ${FILE}"
  exit 1
fi

INPUT="${FILEDIR}/${FILE}"
if [[ ! -f "${INPUT}" ]]; then
  echo "Input file not found on disk: ${INPUT}"
  exit 1
fi

echo "Using input: ${INPUT}"
echo

FAIL=0

for FCL in "${FHICLS[@]}"; do
  NAME="$(basename "${FCL}" .fcl)"
  LOG="test_${NAME}.log"

  echo "==> TEST ${FCL}"

  if [[ ! -f "${FCL}" ]]; then
    echo "MISSING ${FCL}"
    FAIL=1
    echo
    continue
  fi

  rm -f "${LOG}"

  if lar -c "${FCL}" -s "${INPUT}" -n "${NUM_EVENTS}" > "${LOG}" 2>&1; then
    echo "PASS ${FCL}"
  else
    echo "FAIL ${FCL}  (see ${LOG})"
    FAIL=1
  fi

  echo
done

if [[ "${FAIL}" -ne 0 ]]; then
  exit 1
fi

echo "All FHiCL tests passed."
exit 0
