#!/bin/sh
# submit_jobs.sh

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: source submit_grid_jobs.sh <xml_config_file> <stage>"
    exit 1
fi

xml_config_file=$1
stage=$2

tarball_name="StrangenessCode.tar"
tarball_dest="/pnfs/uboone/resilient/users/nlane/NeutralKaon/tarballs/"
scratch_log="/pnfs/uboone/scratch/users/nlane/kaon_dl/v08_00_00_83/nlane_k0signal_overlay_training_nohadrons_reco2_reco2/make_csv/log/jobids.list"
grid_link="https://landscape.fnal.gov/lens/view/job/15660692.0@jobsub05.fnal.gov/"

BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
GREEN="\033[1;32m"
DEFAULT="\033[0m"

echo -e "${BLUE}Starting grid job submission process...${DEFAULT}"

echo -e "${BLUE}Creating tarball of the current code...${DEFAULT}"
cd $MRB_TOP || { echo -e "${RED}Failed to navigate to MRB_TOP${DEFAULT}"; exit 1; }
make_tar_uboone.sh $tarball_name

echo -e "${BLUE}Copying tarball to resilient directory...${DEFAULT}"
cp $tarball_name $tarball_dest || { echo -e "${RED}Failed to copy tarball${DEFAULT}"; exit 1; }
echo -e "${GREEN}Tarball copied successfully.${DEFAULT}"

echo -e "${BLUE}Authenticating with htgettoken...${DEFAULT}"
htgettoken -a htvaultprod.fnal.gov -i uboone || { echo -e "${RED}Authentication failed!${DEFAULT}"; exit 1; }
echo -e "${GREEN}Authentication successful.${DEFAULT}"

echo -e "${BLUE}Submitting jobs with XML configuration: $xml_config_file and stage: $stage...${DEFAULT}"
project.py --xml $xml_config_file --stage $stage --submit

echo -e "${BLUE}Displaying job logs...${DEFAULT}"
if [ -f "$scratch_log" ]; then
    cat "$scratch_log"
else
    echo -e "${RED}Job log not found at: $scratch_log${DEFAULT}"
fi

echo -e "${YELLOW}You can monitor your grid jobs at:${DEFAULT} ${BLUE}$grid_link${DEFAULT}"