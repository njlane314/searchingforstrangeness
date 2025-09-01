#!/bin/bash
set -ex

echo "--- Checking for CVMFS ---"
if [ ! -d /cvmfs ]; then
  echo "Error: /cvmfs is not accessible inside the container." >&2
  echo "You may need to explicitly bind it, e.g., 'apptainer exec --bind /cvmfs ...'." >&2
  exit 1
fi
echo "--- CVMFS found ---"

echo "--- Sourcing minimal UPS setup scripts ---"
# This is the direct alternative to setup_uboone.sh for an AL9/Ubuntu-like environment
# Removed error suppression to debug script failure
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
source /cvmfs/uboone.opensciencegrid.org/products/setup
source /cvmfs/larsoft.opensciencegrid.org/products/setup
echo "--- Minimal UPS setup complete ---"

echo "--- Setting up hdf5 (suppressing harmless warnings) ---"
# Redirect stderr to /dev/null for this command to hide known, non-fatal errors
setup hdf5 v1_12_2a -q e20:prof 2>/dev/null
echo "--- hdf5 setup complete ---"


echo "--- Locating and sourcing ROOT ---"
if command -v root-config >/dev/null 2>&1; then
  source "$(root-config --prefix)/bin/thisroot.sh"
else
  source "/usr/local/root/bin/thisroot.sh"
fi
echo "--- ROOT setup complete ---"

INPUT_FILE="$1"
OUTPUT_FILE="$2"
WEIGHTS_FILE="$3"
TREE_NAME="$4"
BRANCH_NAME="$5"

echo "--- Input arguments received ---"
echo "Input file: ${INPUT_FILE}"
echo "Output file: ${OUTPUT_FILE}"
echo "Weights file: ${WEIGHTS_FILE}"
echo "Tree name: ${TREE_NAME}"
echo "Branch name: ${BRANCH_NAME}"

unset PYTHONHOME
unset PYTHONPATH

echo "--- Setting up Python from CVMFS ---"
PY_SETUP=$(find /cvmfs/uboone.opensciencegrid.org/products/python /cvmfs/larsoft.opensciencegrid.org/products/python -name "setup.sh" 2>/dev/null | head -n 1)
if [ -n "$PY_SETUP" ]; then
    source "$PY_SETUP"
    echo "--- Sourced Python from $PY_SETUP ---"
else
    echo "--- Could not find Python setup script in CVMFS ---"
fi

echo "--- Starting Python inference script ---"
python3 run_inference.py \
  --input "$INPUT_FILE" \
  --output "$OUTPUT_FILE" \
  --weights "$WEIGHTS_FILE" \
  --tree "$TREE_NAME" \
  --branch "$BRANCH_NAME"
echo "--- Python script finished ---"

echo "--- Verifying output file creation ---"
if [ ! -f "$OUTPUT_FILE" ]; then
  echo "Error: Python script did not create the expected output file: ${OUTPUT_FILE}" >&2
  exit 1
fi
echo "--- Script completed successfully, output file found. ---"
