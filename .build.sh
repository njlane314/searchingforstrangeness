#!/bin/bash

export WRK_DIR=$(pwd)
cd $MRB_TOP

mrbsetenv
mrb i -j10 

cd $WRK_DIR