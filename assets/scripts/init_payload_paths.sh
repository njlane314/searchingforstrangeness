#!/bin/bash
set -euo pipefail

# FHiCLs are in your local build tarball:
FCL_DIR="${MRB_INSTALL:-$PWD}/fcl"
[ -d "$FCL_DIR" ] && export FHICL_FILE_PATH="$FCL_DIR:${FHICL_FILE_PATH:-}"

# Assets (weights/models/scripts) come from the RCDS tarball:
ASSETS="${INPUT_TAR_DIR_LOCAL:-}"
# Fallback: if you ever tuck assets into your local build:
[ -z "$ASSETS" ] && ASSETS="${MRB_INSTALL:-$PWD}"

# Discovery for modules that use cet::search_path("FW_SEARCH_PATH")
export FW_SEARCH_PATH="${ASSETS}:${ASSETS}/weights:${ASSETS}/models:${FW_SEARCH_PATH:-}"

# Convenience vars (handy in FHiCL with %(ENV{...}))
export WEIGHTS_DIR="${WEIGHTS_DIR:-${ASSETS}/weights}"
export MODELS_DIR="${MODELS_DIR:-${ASSETS}/models}"
export SCRIPTS_DIR="${SCRIPTS_DIR:-${ASSETS}/scripts}"

# Optional helpers in the assets tarball
[ -d "${ASSETS}/scripts" ] && export PATH="${ASSETS}/scripts:${PATH}"
[ -d "${ASSETS}/python"  ] && export PYTHONPATH="${ASSETS}/python:${PYTHONPATH:-}"

echo "FHICL_FILE_PATH=${FHICL_FILE_PATH}"
echo "FW_SEARCH_PATH=${FW_SEARCH_PATH}"
