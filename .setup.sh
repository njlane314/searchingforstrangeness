#!/bin/bash
set -euo pipefail

set -x

UBOONECODE_VERSION="${UBOONECODE_VERSION:-v08_00_00_82}"
PYTHON_VERSION="${PYTHON_VERSION:-v2_7_14b}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --uboonecode-version) UBOONECODE_VERSION="${2:?--uboonecode-version requires an argument}"; shift 2 ;;
    --python-version)     PYTHON_VERSION="${2:?--python-version requires an argument}"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

unsetup_all || true
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode "${UBOONECODE_VERSION}" -q "e17:prof"
setup sam_web_client

if ! command -v mrb >/dev/null 2>&1; then
  setup mrb || { echo "Failed to setup mrb" >&2; exit 1; }
fi

shopt -s nullglob
matches=(../../../../localProducts_*/setup)
shopt -u nullglob
if (( ${#matches[@]} == 0 )); then
  echo "No localProducts_* setup file found under ../../../../" >&2
  exit 1
elif (( ${#matches[@]} > 1 )); then
  echo "Multiple localProducts setups found: ${matches[*]}" >&2
  echo "Please specify the correct one." >&2
  exit 1
else
  source "${matches[0]}"
fi

if ! mrbslp; then
  echo "mrbslp failed" >&2
  exit 1
fi

if ! command -v htgettoken >/dev/null 2>&1; then
  echo "htgettoken not found in PATH" >&2
  exit 1
fi
if ! htgettoken -a htvaultprod.fnal.gov -i uboone; then
  echo "htgettoken failed" >&2
  exit 1
fi

if ! command -v jobsub_submit >/dev/null 2>&1; then
  echo "jobsub_submit not found in PATH (continuing)" >&2
fi

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}/cvmfs/larsoft.opensciencegrid.org/products/python/${PYTHON_VERSION}/Linux64bit+3.10-2.17/lib"

if ! command -v kx509 >/dev/null 2>&1; then
  echo "kx509 not found in PATH" >&2
  exit 1
fi
if ! kx509; then
  echo "kx509 failed" >&2
  exit 1
fi
if ! voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis; then
  echo "voms-proxy-init failed" >&2
  exit 1
fi

echo "Environment setup complete."

