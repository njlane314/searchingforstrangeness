#!/bin/bash
set -e

# Locate and source the ROOT setup script. Prefer a standard installation
# discovered via `root-config` but fall back to the historical location
# if that fails.  Failing to source ROOT should terminate with a clear
# error message so the calling art job reports the problem.
if command -v root-config >/dev/null 2>&1; then
  ROOT_PREFIX="$(root-config --prefix)"
  ROOT_SETUP="${ROOT_PREFIX}/bin/thisroot.sh"
elif [ -f /usr/local/root/bin/thisroot.sh ]; then
  ROOT_SETUP="/usr/local/root/bin/thisroot.sh"
else
  echo "Error: unable to locate ROOT setup script" >&2
  exit 1
fi
source "$ROOT_SETUP"
INPUT_FILE="$1"
OUTPUT_FILE="$2"
WEIGHTS_FILE="$3"
TREE_NAME="$4"
BRANCH_NAME="$5"
unset PYTHONHOME
unset PYTHONPATH
# Attempt to source a Python environment from CVMFS so that required
# libraries like h5py are available.
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
  echo "Warning: No CVMFS Python setup script found; using system Python." >&2
fi
if ! python3 -c "import h5py" >/dev/null 2>&1; then
  echo "Error: Python module 'h5py' is required but not installed." >&2
  exit 1
fi
# Check for additional Python dependencies required by run_inference.py.
if ! python3 -c "import MinkowskiEngine" >/dev/null 2>&1; then
  echo "Error: Python module 'MinkowskiEngine' is required but not installed." >&2
  exit 1
fi
if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file '$INPUT_FILE' is missing or not provided." >&2
  exit 1
fi
if [ -z "$WEIGHTS_FILE" ] || [ ! -f "$WEIGHTS_FILE" ]; then
  echo "Error: Weights file '$WEIGHTS_FILE' is missing or not provided." >&2
  exit 1
fi
echo "Input file: ${INPUT_FILE}"
echo "Output file: ${OUTPUT_FILE}"
echo "Weights file: ${WEIGHTS_FILE}"
echo "Tree name: ${TREE_NAME}"
echo "Branch name: ${BRANCH_NAME}"
python3 run_inference.py \
  --input "$INPUT_FILE" \
  --output "$OUTPUT_FILE" \
  --weights "$WEIGHTS_FILE" \
  --tree "$TREE_NAME" \
  --branch "$BRANCH_NAME"

# Ensure the inference produced an output file; report the error to stderr so
# that the art framework captures the reason for a non-zero exit status.
if [ ! -f "$OUTPUT_FILE" ]; then
  echo "Error: Python script did not create the expected output file: ${OUTPUT_FILE}" >&2
  exit 1
fi

