#!/bin/bash
set -euo pipefail


FCL_DIR="${MRB_INSTALL:-$PWD}/fcl"
[ -d "$FCL_DIR" ] && export FHICL_FILE_PATH="$FCL_DIR:${FHICL_FILE_PATH:-}"


ASSETS="${INPUT_TAR_DIR_LOCAL:-}"

[ -z "$ASSETS" ] && ASSETS="${MRB_INSTALL:-$PWD}"


export FW_SEARCH_PATH="${ASSETS}:${ASSETS}/weights:${ASSETS}/models:${FW_SEARCH_PATH:-}"


export WEIGHTS_DIR="${WEIGHTS_DIR:-${ASSETS}/weights}"
export MODELS_DIR="${MODELS_DIR:-${ASSETS}/models}"
export SCRIPTS_DIR="${SCRIPTS_DIR:-${ASSETS}/scripts}"


[ -d "${ASSETS}/scripts" ] && export PATH="${ASSETS}/scripts:${PATH}"
[ -d "${ASSETS}/python"  ] && export PYTHONPATH="${ASSETS}/python:${PYTHONPATH:-}"

echo "FHICL_FILE_PATH=${FHICL_FILE_PATH}"
echo "FW_SEARCH_PATH=${FW_SEARCH_PATH}"

