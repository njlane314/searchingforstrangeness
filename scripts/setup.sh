unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup larsoft v10_04_07 -q e26:prof
setup libtorch v2_1_1a -q e26
setup mrb

source /exp/uboone/app/users/nlane/production/KaonShortProduction05/localProducts_larsoft_v10_04_07_e26_prof/setup
mrbslp

alias cdt="cd ${MRB_TOP}"

export SEARCH_TOP="/exp/uboone/app/users/nlane/production/KaonShortProduction05/srcs/ubana/ubana/searchingforstrangeness"
alias cds="cd ${SEARCH_TOP}"

export TORCH_DIR=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v2_1_1a/Linux64bit+3.10-2.17-e26/lib/python3.8/site-packages/torch/share/cmake/Torch
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v3_8_10/Linux64bit+3.10-2.17/lib
