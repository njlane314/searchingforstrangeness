#!/bin/sh

set -e

BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
GREEN="\033[1;32m"
DEFAULT="\033[0m"

USER=$(echo $USER)

if [ "$#" -ne 2 ] && [ "$#" -ne 0 ]; then
    echo "Usage: source check_jobs.sh [xml_config_file] [stage]"
    return 1
fi

if [ "$#" -eq 2 ]; then
    xml_config_file=$1
    stage=$2

    echo -e "${BLUE}-- Checking job status for XML: $xml_config_file and Stage: $stage...${DEFAULT}"

    log_dir=$(grep -oP '(?<=<logdir>).*?(?=</logdir>)' "$xml_config_file")
    release=$(grep -oP '(?<=<!ENTITY release ").*?(?=">)' "$xml_config_file")
    project_name=$(grep -oP '(?<=<!ENTITY name ").*?(?=">)' "$xml_config_file")

    log_dir=$(echo "$log_dir" | sed "s/&release;/$release/" | sed "s/&name;/$project_name/")
    log_file="${log_dir}/jobids.list"

    if [ -f "$log_file" ]; then
        echo -e "${BLUE}-- Reading job ID from log file: $log_file...${DEFAULT}"
        job_id=$(cat "$log_file" | head -n 1)

        if [ -n "$job_id" ]; then
            echo -e "${BLUE}-- Found job ID: $job_id${DEFAULT}"
            echo -e "${BLUE}-- Link: https://landscape.fnal.gov/lens/view/job/${job_id}${DEFAULT}"
        else
            echo -e "${RED}-- No job ID found in $log_file.${DEFAULT}"
        fi
    else
        echo -e "${RED}-- Log file not found: $log_file.${DEFAULT}"
    fi

    echo -e "${GREEN}-- Specific job check complete!${DEFAULT}"
    return 0
fi

if [ "$#" -eq 0 ]; then
    echo -e "${BLUE}-- Checking all jobs for user '${USER}'...${DEFAULT}"
    job_queue=$(jobsub_q $USER)

    echo -e "${BLUE}-- Job queue output:${DEFAULT}"
    echo "$job_queue"

    echo -e "${YELLOW}-- Generating detailed monitoring links for all jobs...${DEFAULT}"
    echo "$job_queue" | awk '/@/ {print $1}' | while read -r jobid; do
        echo -e "${BLUE}-- Link: https://landscape.fnal.gov/lens/view/job/${jobid}${DEFAULT}"
    done

    echo -e "${GREEN}-- General job check complete!${DEFAULT}"
    return 0
fi

echo -e "${GREEN}-- Job status check complete!${DEFAULT}"