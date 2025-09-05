#!/usr/bin/env bash

TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"
BUILD_ROOT="${BUILD_ROOT:-${MRB_INSTALL:-${SRT_PRIVATE_CONTEXT:-${SCRIPT_DIR}}}}"

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

LAR_TAR_LOCAL="${TMP_DIR}/strangeness.tar"
LAR_TAR="${TAR_DIR}/strangeness.tar"

tar -C "${BUILD_ROOT}" --exclude='.git' --exclude='tmp' --exclude='*.root' -czf "${LAR_TAR_LOCAL}" .

mkdir -p "${TAR_DIR}"
mv -f "${LAR_TAR_LOCAL}" "${LAR_TAR}"

echo "LArSoft tarball:  ${LAR_TAR}"
