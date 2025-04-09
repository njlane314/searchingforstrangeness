#!/bin/sh
# submit_jobs.sh
#project.py --checkana --xml scripts/xmls/make_trainingsamples.xml --stage analyse

set -e

BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
GREEN="\033[1;32m"
DEFAULT="\033[0m"

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "-- Usage: source submit_grid_jobs.sh <xml_config_file> <stage> [skip-tar]"
    return 1
fi

xml_config_file=$1
stage=$2
skip_tarball=${3:-false}

if [ -z "$TARBALL_NAME" ] || [ -z "$TARBALL_DEST" ]; then
    echo -e "${RED}-- Error: TARBALL_NAME or TARBALL_DEST is not set. Please source the setup script.${DEFAULT}"
    return 1
fi

echo -e "${BLUE}-- Starting grid job submission process...${DEFAULT}"

if [ "$skip_tarball" != "skip-tar" ]; then
    echo -e "${BLUE}-- Creating tarball of the current code...${DEFAULT}"
    cd $MRB_TOP || { echo -e "${RED}-- Failed to navigate to MRB_TOP${DEFAULT}"; return 1; }
    make_tar_uboone.sh $TARBALL_NAME

    echo -e "${BLUE}-- Copying tarball to resilient directory...${DEFAULT}"
    cp -f $TARBALL_NAME $TARBALL_DEST || { echo -e "${RED}-- Failed to copy tarball${DEFAULT}"; return 1; }
    echo -e "${GREEN}-- Tarball copied successfully to $TARBALL_DEST.${DEFAULT}"
else
    echo -e "${YELLOW}-- Skipping tarball creation and copy step as requested.${DEFAULT}"
fi

echo -e "${BLUE}-- Authenticating with htgettoken...${DEFAULT}"
htgettoken -a htvaultprod.fnal.gov -i uboone || { echo -e "${RED}-- Authentication failed!${DEFAULT}"; return 1; }
echo -e "${GREEN}-- Authentication successful.${DEFAULT}"

echo -e "${BLUE}-- Submitting jobs with XML configuration: $xml_config_file and stage: $stage...${DEFAULT}"
cd ${SEARCH_TOP}
if ! project.py --xml $xml_config_file --stage $stage --submit; then
    echo -e "${YELLOW}-- Submission failed. Cleaning stage and retrying...${DEFAULT}"
    
    if project.py --xml $xml_config_file --stage $stage --clean; then
        echo -e "${BLUE}-- Stage cleaned successfully. Retrying submission...${DEFAULT}"
        
        if ! project.py --xml $xml_config_file --stage $stage --submit; then
            echo -e "${RED}-- Submission failed again after cleaning. Please check for issues.${DEFAULT}"
            return 1
        else
            echo -e "${GREEN}-- Submission succeeded after cleaning.${DEFAULT}"
        fi
    else
        echo -e "${RED}-- Failed to clean the stage. Please investigate manually.${DEFAULT}"
        return 1
    fi
else
    echo -e "${GREEN}-- Submission succeeded on the first attempt.${DEFAULT}"
fi
