#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env/setenv.sh"
source "${SCRIPT_DIR}/env/configure.sh"
unset PYTHONHOME
unset PYTHONPATH
# Try to source a Python setup script from CVMFS to match the runtime
# environment used by the inference wrapper.
PY_SETUP=""
for pyroot in \
  /cvmfs/uboone.opensciencegrid.org/products/python \
  /cvmfs/larsoft.opensciencegrid.org/products/python; do
  if [ -d "$pyroot" ]; then
    PY_SETUP=$(find "$pyroot" -maxdepth 2 -name setup.sh 2>/dev/null | head -n 1)
    if [ -n "$PY_SETUP" ]; then
      source "$PY_SETUP"
      break
    fi
  fi
done
if [ -z "$PY_SETUP" ]; then
  echo "Warning: No CVMFS Python setup script found; using system Python."
fi
echo "--- Starting container test script (test_container_env.sh) ---"
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -S) INPUT_LIST="$2"; shift; shift;;
    --outdir) JOB_OUTDIR="$2"; shift; shift;;
    *) shift;;
  esac
done
echo "This is a TEST job."
echo "Input file list is: $INPUT_LIST"
echo "Output directory is: $JOB_OUTDIR"
echo "========================================"
echo
echo "--- 1. Checking OS release information ---"
cat /etc/os-release
echo
echo "--- 2. Checking Python version and location ---"
which python3
python3 --version
echo
echo "--- 3. Checking for PyTorch ---"
python3 -c "import torch; print('PyTorch version:', torch.__version__)"
echo
echo "--- 4. Checking for uproot ---"
python3 -c "import uproot; print('uproot version:', uproot.__version__)"
echo
echo "--- 5. Checking for h5py ---"
python3 -c "import h5py; print('h5py version:', h5py.__version__)" || echo "h5py not available"
echo "--- Test script finished ---"
exit 0

