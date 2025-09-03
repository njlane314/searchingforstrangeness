#!/usr/bin/env bash
#
# .tar.sh — build and publish tarballs for local larsoft and analysis assets
#
# This script creates two tarballs and places them in a common tarball
# directory:
#   1. LArSoft build (strangeness.tar)
#   2. Analysis assets (strangeness_assets.tar.gz)
#
# Environment variables used for customisation:
#   TAR_DIR      - destination directory for the tarballs
#   ASSETS_ROOT  - location of the assets tree
#
# After running, the following variables are exported for convenience:
#   STRANGENESS_TAR         - full path to the LArSoft tarball
#   STRANGENESS_ASSETS_TAR  - full path to the assets tarball

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Destination for tarballs; default to user's scratch area
TAR_DIR="${TAR_DIR:-/pnfs/uboone/scratch/users/${USER}/tarballs}"
# Location of the assets to package
ASSETS_ROOT="${ASSETS_ROOT:-${REPO_ROOT}/assets}"

LARSOFT_TAR="${TAR_DIR}/strangeness.tar"
ASSETS_TAR="${TAR_DIR}/strangeness_assets.tar.gz"

mkdir -p "${TAR_DIR}"

pushd "${REPO_ROOT}" >/dev/null

# Build the larsoft tarball
"${SCRIPT_DIR}/.tar_uboone.sh" strangeness.tar
mv -f strangeness.tar "${LARSOFT_TAR}"

# Build the assets tarball
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

