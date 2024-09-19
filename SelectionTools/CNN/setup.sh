#!/bin/bash

source /cvmfs/larsoft.opensciencegrid.org/setup_larsoft.sh

setup cmake v3_24_1 -q Linux64bit+3.10-2.17
setup gcc v12_1_0 -q Linux64bit+3.10-2.17
setup libtorch v2_1_1b -q Linux64bit+3.10-2.17:e26
setup caffe v1_0k -q Linux64bit+3.10-2.17:e17:prof
setup tbb v2021_9_0 -q Linux64bit+3.10-2.17:e26
setup protobuf v3_21_12 -q Linux64bit+3.10-2.17:e26
setup root v6_28_12 -q Linux64bit+3.10-2.17:e26:p3915:prof