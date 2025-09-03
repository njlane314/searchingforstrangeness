#!/usr/bin/env bash
















set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"


TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"

ASSETS_ROOT="${ASSETS_ROOT:-${REPO_ROOT}/assets}"

LARSOFT_TAR="${TAR_DIR}/strangeness.tar"
ASSETS_TAR="${TAR_DIR}/strangeness_assets.tar.gz"

mkdir -p "${TAR_DIR}"

pushd "${REPO_ROOT}" >/dev/null


"${SCRIPT_DIR}/.tar_uboone.sh" strangeness.tar
mv -f strangeness.tar "${LARSOFT_TAR}"


"${SCRIPT_DIR}/.pack_assets.sh" -r "${ASSETS_ROOT}" -o strangeness_assets.tar.gz
mv -f strangeness_assets.tar.gz "${ASSETS_TAR}"

popd >/dev/null

export STRANGENESS_TAR="${LARSOFT_TAR}"
export STRANGENESS_ASSETS_TAR="${ASSETS_TAR}"

echo "LArSoft tarball:  ${STRANGENESS_TAR}"
echo "Assets tarball:  ${STRANGENESS_ASSETS_TAR}"
echo "Environment variables set:"
echo "  export STRANGENESS_TAR=${STRANGENESS_TAR}"
echo "  export STRANGENESS_ASSETS_TAR=${STRANGENESS_ASSETS_TAR}"

