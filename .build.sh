#!/bin/bash

set -euo pipefail

# Ensure required environment variables and commands are available
if [[ -z "${MRB_TOP:-}" ]]; then
  echo "MRB_TOP is not set; run .setup.sh before building." >&2
  exit 1
fi

if ! command -v mrbsetenv >/dev/null 2>&1; then
  echo "mrbsetenv not found; ensure the mrb environment is configured." >&2
  exit 1
fi

WRK_DIR=$(pwd)
cd "$MRB_TOP"

mrbsetenv
mrb i -j8

cd "$WRK_DIR"

# Package the LArSoft local installation and assets if the build succeeds
"$WRK_DIR/scripts/tar.sh"
