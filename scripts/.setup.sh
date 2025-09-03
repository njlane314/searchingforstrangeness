#!/bin/bash

# Load the UPS setup to make `setup` and `unsetup_all` available.
if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

# Ensure MRB is available before attempting to use it.
if ! command -v mrb >/dev/null 2>&1; then
  setup mrb
fi

# Reset any existing UPS products to avoid environment conflicts.
unsetup_all || true

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_82 -q e17:prof
#setup libtorch v1_0_1 -q Linux64bit+3.10-2.17:e17:prof
setup hdf5 v1_10_5 -q e17

setup sam_web_client

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbslp

htgettoken -a htvaultprod.fnal.gov -i uboone

which jobsub_submit

export TORCH_DIR_BASE=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1/Linux64bit+3.10-2.17-e17-prof/lib/python2.7
export TORCH_DIR="$TORCH_DIR_BASE/site-packages/torch/share/cmake/Torch"
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib"

kx509
voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis

