source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_82 -q "e17:prof"

unsetup mrb
setup mrb -o

source /exp/uboone/app/users/nlane/production/KaonShortProduction01/localProducts_larsoft_v08_05_00_24_e17_prof/setup
mrbslp

alias cdt="cd ${MRB_TOP}"