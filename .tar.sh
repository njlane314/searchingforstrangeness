#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"
ASSETS_ROOT="${ASSETS_ROOT:-${SCRIPT_DIR}/assets}"
BUILD_ROOT="${BUILD_ROOT:-${MRB_INSTALL:-${SRT_PRIVATE_CONTEXT:-${SCRIPT_DIR}}}}"

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

LAR_TAR_LOCAL="${TMP_DIR}/strangeness.tar"
ASSETS_TAR_LOCAL="${TMP_DIR}/strangeness_assets.tar.gz"
LAR_TAR="${TAR_DIR}/strangeness.tar"
ASSETS_TAR="${TAR_DIR}/strangeness_assets.tar.gz"

tar -C "${BUILD_ROOT}" --exclude='.git' --exclude='tmp' --exclude='*.root' -czf "${LAR_TAR_LOCAL}" .
tar -C "${ASSETS_ROOT}" --exclude='.git' --exclude='__pycache__' --exclude='*.pyc' -czf "${ASSETS_TAR_LOCAL}" .

mkdir -p "${TAR_DIR}"
mv -f "${LAR_TAR_LOCAL}" "${LAR_TAR}"
mv -f "${ASSETS_TAR_LOCAL}" "${ASSETS_TAR}"

echo "LArSoft tarball:  ${LAR_TAR}"
echo "Assets tarball:  ${ASSETS_TAR}"
