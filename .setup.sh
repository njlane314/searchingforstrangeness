#!/bin/bash

unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup larsoft v10_04_07 -q e26:prof
setup uboonecode v10_04_07_04 -q e26:prof
#setup libtorch v2_1_1a -q e26
setup sam_web_client
setup mrb

source /exp/uboone/app/users/nlane/production/strangeness_mcc10/localProducts_*/setup
mrbslp

htgettoken -a htvaultprod.fnal.gov -i uboone

