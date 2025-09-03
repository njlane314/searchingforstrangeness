#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"
ASSETS_ROOT="${ASSETS_ROOT:-${SCRIPT_DIR}/assets}"
BUILD_ROOT="${BUILD_ROOT:-${MRB_INSTALL:-${SRT_PRIVATE_CONTEXT:-${SCRIPT_DIR}}}}"

mkdir -p "${TAR_DIR}"

LAR_TAR="${TAR_DIR}/strangeness.tar"
ASSETS_TAR="${TAR_DIR}/strangeness_assets.tar.gz"

if [[ -d "${BUILD_ROOT}" ]]; then
  tar -C "${BUILD_ROOT}" -czf "${LAR_TAR}" . \
    --exclude='.git' --exclude='tmp' --exclude='*.root'
else
  echo "WARNING: build directory '${BUILD_ROOT}' not found" >&2
fi

if [[ -d "${ASSETS_ROOT}" ]]; then
  tar -C "${ASSETS_ROOT}" -czf "${ASSETS_TAR}" . \
    --exclude='.git' --exclude='__pycache__' --exclude='*.pyc'
else
  echo "WARNING: assets directory '${ASSETS_ROOT}' not found" >&2
fi

echo "LArSoft tarball:  ${LAR_TAR}"
echo "Assets tarball:  ${ASSETS_TAR}"
