#!/bin/bash
export WRK_DIR=$(pwd)
cd $MRB_TOP

mrbsetenv
mrb i -j8

cd $WRK_DIR
source $WRK_DIR/.tar.sh

source assets/script/init_payload_paths.sh
