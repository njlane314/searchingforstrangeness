#!/bin/bash
set -euo pipefail

# Default software versions. Override via environment variables or
# command-line arguments (e.g., --uboonecode-version vX_Y_Z).
UBOONECODE_VERSION="${UBOONECODE_VERSION:-v08_00_00_82}"
HDF5_VERSION="${HDF5_VERSION:-v1_10_5}"
LIBTORCH_VERSION="${LIBTORCH_VERSION:-v1_0_1}"
PYTHON_VERSION="${PYTHON_VERSION:-v2_7_14b}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --uboonecode-version)
      UBOONECODE_VERSION="${2:?--uboonecode-version requires an argument}"
      shift 2
      ;;
    --hdf5-version)
      HDF5_VERSION="${2:?--hdf5-version requires an argument}"
      shift 2
      ;;
    --libtorch-version)
      LIBTORCH_VERSION="${2:?--libtorch-version requires an argument}"
      shift 2
      ;;
    --python-version)
      PYTHON_VERSION="${2:?--python-version requires an argument}"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Load the UPS setup to make `setup` and `unsetup_all` available.
if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

# Ensure MRB is available before attempting to use it.
if ! command -v mrb >/dev/null 2>&1; then
  setup mrb
fi

# Reset any existing UPS products to avoid environment conflicts.
unsetup_all || true

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode "$UBOONECODE_VERSION" -q e17:prof
#setup libtorch "$LIBTORCH_VERSION" -q Linux64bit+3.10-2.17:e17:prof
setup hdf5 "$HDF5_VERSION" -q e17

setup sam_web_client

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
if ! mrbslp; then
  echo "mrbslp failed" >&2
  exit 1
fi

if ! htgettoken -a htvaultprod.fnal.gov -i uboone; then
  echo "htgettoken failed" >&2
  exit 1
fi

which jobsub_submit

export TORCH_DIR_BASE="/cvmfs/uboone.opensciencegrid.org/products/libtorch/${LIBTORCH_VERSION}/Linux64bit+3.10-2.17-e17-prof/lib/python2.7"
export TORCH_DIR="$TORCH_DIR_BASE/site-packages/torch/share/cmake/Torch"
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/${PYTHON_VERSION}/Linux64bit+3.10-2.17/lib"

if ! kx509; then
  echo "kx509 failed" >&2
  exit 1
fi
voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis

