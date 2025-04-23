#!/bin/bash

unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup larsoft v10_04_07 -q e26:prof
setup uboonecode v10_04_07_04 -q e26:prof
#setup uboonecode v08_00_00_88 -q e17:prof
setup sam_web_client
setup mrb

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbslp

htgettoken -a htvaultprod.fnal.gov -i uboone

