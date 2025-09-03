#!/bin/bash
set -e

export WRK_DIR=$(pwd)
cd "$MRB_TOP"

mrbsetenv
mrb i -j8

cd "$WRK_DIR"
"$WRK_DIR/tar_artifacts.sh"
