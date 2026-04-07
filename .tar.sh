#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"
BUILD_ROOT="${BUILD_ROOT:-${MRB_INSTALL:-${SRT_PRIVATE_CONTEXT:-${SCRIPT_DIR}}}}"
REINT_DATA_BUILD_DIR="${BUILD_ROOT}/job/reinteractions/data"
REINT_DATA_SOURCE_DIR="${SCRIPT_DIR}/job/reinteractions/data"

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

LAR_TAR_LOCAL="${TMP_DIR}/strangeness.tar"
LAR_TAR="${TAR_DIR}/strangeness.tar"
STAGE_DIR="${TMP_DIR}/stage"
REINT_DATA_DIR=""

if [ -d "${REINT_DATA_BUILD_DIR}" ]; then
  REINT_DATA_DIR="${REINT_DATA_BUILD_DIR}"
elif [ -d "${REINT_DATA_SOURCE_DIR}" ]; then
  REINT_DATA_DIR="${REINT_DATA_SOURCE_DIR}"
fi

mkdir -p "${STAGE_DIR}"
tar -C "${BUILD_ROOT}" --exclude='.git' --exclude='tmp' --exclude='*.root' -cf - . | tar -C "${STAGE_DIR}" -xf -

if [ -n "${REINT_DATA_DIR}" ]; then
  mkdir -p "${STAGE_DIR}/job/reinteractions"
  cp -R "${REINT_DATA_DIR}" "${STAGE_DIR}/job/reinteractions/"
fi

tar -C "${STAGE_DIR}" -czf "${LAR_TAR_LOCAL}" .

mkdir -p "${TAR_DIR}"
mv -f "${LAR_TAR_LOCAL}" "${LAR_TAR}"

echo "LArSoft tarball:  ${LAR_TAR}"
