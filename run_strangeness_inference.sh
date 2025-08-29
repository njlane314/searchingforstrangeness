#!/bin/bash
set -e
source /usr/local/root/bin/thisroot.sh
INPUT_FILE="$1"
OUTPUT_FILE="$2"
WEIGHTS_FILE="$3"
TREE_NAME="$4"
BRANCH_NAME="$5"
unset PYTHONHOME
unset PYTHONPATH
if ! python3 -c "import h5py" >/dev/null 2>&1; then
  echo "Error: Python module 'h5py' is required but not installed."
  exit 1
fi
if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file '$INPUT_FILE' is missing or not provided."
  exit 1
fi
if [ -z "$WEIGHTS_FILE" ] || [ ! -f "$WEIGHTS_FILE" ]; then
  echo "Error: Weights file '$WEIGHTS_FILE' is missing or not provided."
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
if [ ! -f "$OUTPUT_FILE" ]; then
  echo "Error: Python script did not create the expected output file: ${OUTPUT_FILE}"
  exit 1
fi

