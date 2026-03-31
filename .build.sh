#!/bin/bash
export WRK_DIR=$(pwd)
cd $MRB_TOP

mrbsetenv
# ROOT dictionary checks resolve headers through ROOT_INCLUDE_PATH rather than
# the compiler include flags used during the build.
export ROOT_INCLUDE_PATH="$WRK_DIR${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
mrb i -j8

cd $WRK_DIR
source $WRK_DIR/.tar.sh
