#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env/setenv.sh"
source "${SCRIPT_DIR}/env/configure.sh"
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
echo "--- Test script finished ---"
exit 0

