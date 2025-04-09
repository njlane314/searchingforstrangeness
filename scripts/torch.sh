#!/bin/bash

setup libtorch v1_0_1 -q Linux64bit+3.10-2.17:e17:prof

export TORCH_DIR=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1/Linux64bit+3.10-2.17-e17-prof/lib/python2.7/site-packages/torch/share/cmake/Torch
#export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib
