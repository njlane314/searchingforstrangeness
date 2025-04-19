unsetup_all
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup larsoft v10_04_07 -q e26:prof
#setup libtorch v2_1_1a -q e26
setup mrb

source /exp/uboone/app/users/nlane/production/strangeness_mcc10/localProducts_*/setup
mrbslp
