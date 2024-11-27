unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_82 -q "e17:prof"

#ups list -aK+ libtorch
#setup libtorch v1_0_1_mkldnn -q Linux64bit+3.10-2.17:e17:prof
#setup libtorch v1_0_1_nonuma -q Linux64bit+3.10-2.17:e17:prof
setup libtorch v1_0_1 -q Linux64bit+3.10-2.17:e17:prof

unsetup mrb
setup mrb -o

source /exp/uboone/app/users/nlane/production/KaonShortProduction04/localProducts_larsoft_v08_05_00_24_e17_prof/setup
mrbslp

alias cdt="cd ${MRB_TOP}"

export SEARCH_TOP="/exp/uboone/app/users/nlane/production/KaonShortProduction04/srcs/ubana/ubana/searchingforstrangeness"
alias cds="cd ${SEARCH_TOP}"

#export TORCH_DIR=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1_nonuma/Linux64bit+3.10-2.17-e17-prof/share/cmake/Torch/TorchConfig.cmake
export TORCH_DIR=/cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1/Linux64bit+3.10-2.17-e17-prof/lib/python2.7/site-packages/torch/share/cmake/Torch
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$TORCH_DIR"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib

case "$USER" in
    nlane)
        export TARBALL_NAME="StrangenessCode.tar"
        export TARBALL_DEST="/pnfs/uboone/resilient/users/nlane/NeutralKaon/tarballs/"
        ;;
esac

git pull

source scripts/train.sh