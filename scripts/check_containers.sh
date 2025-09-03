#!/usr/bin/env bash
set -euo pipefail

CONTAINER_PATH="/cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod"

echo "Launching container: ${CONTAINER_PATH}"

apptainer shell \
  -B /cvmfs \
  -B /exp/uboone \
  -B /pnfs/uboone \
  -B /run/user \
  -s /bin/bash \
  --env UPS_OVERRIDE="-H Linux64bit+3.10-2.17" \
  "$CONTAINER_PATH" << EOF

echo "Executing Python version checks inside the container."

python3 -c "
import sys
try:
    import torch
    import MinkowskiEngine as ME
    print(f'Python Version: {sys.version.split()[0]}')
    print(f'PyTorch Version: {torch.__version__}')
    print(f'MinkowskiEngine Version: {ME.__version__}')
except ImportError as e:
    print(f'Error: A required package is not installed. {e}', file=sys.stderr)
    exit(1)
"

EOF

echo "Script execution complete."