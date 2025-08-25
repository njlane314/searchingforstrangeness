#!/bin/bash
set -euo pipefail

# Copy required inference files to a PNFS directory so that jobsub can stage them
# to worker nodes with -f options. The destination directory can be provided as
# an argument; otherwise the default resilient path is used.

DEST=${1:-/pnfs/uboone/resilient/users/nlane/misc}

mkdir -p "$DEST"

FILES=(run_strangeness_inference.sh run_inference.py binary_classifier_resnet34.pth badchannels.txt)
for f in "${FILES[@]}"; do
  cp "$f" "$DEST/"
 done

ls "$DEST"
