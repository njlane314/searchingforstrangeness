#!/bin/bash

unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
#setup larsoft v10_04_07 -q e26:prof
#setup uboonecode v10_04_07_04 -q e26:prof
setup uboonecode v08_00_00_82 -q e17:prof
setup libtorch v1_0_1 -q Linux64bit+3.10-2.17:e17:prof

setup sam_web_client
#setup mrb

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbslp

htgettoken -a htvaultprod.fnal.gov -i uboone

which jobsub_submit

export TORCH_DIR=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1/Linux64bit+3.10-2.17-e17-prof/lib/python2.7/site-packages/torch/share/cmake/Torch
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/l
