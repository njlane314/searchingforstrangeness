#!/bin/bash

set -e

export WRK_DIR=$(pwd)
cd "$MRB_TOP"

mrbsetenv
mrb i -j8

cd "$WRK_DIR"

# Package the LArSoft local installation and assets if the build succeeds
"$WRK_DIR/scripts/tar.sh"
