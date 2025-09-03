if ! command -v setup >/dev/null 2>&1; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup
fi

if ! command -v mrb >/dev/null 2>&1; then
  setup mrb
fi

unsetup_all || true

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_82 -q e17:prof

setup sam_web_client

export WRK_DIR=$(pwd)
source ../../../../localProducts_*/setup
mrbslp

htgettoken -a htvaultprod.fnal.gov -i uboone

which jobsub_submit

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib"

kx509
voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis
