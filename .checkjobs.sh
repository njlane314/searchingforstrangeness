#!/bin/sh

set -e

get_xml_entity() {
    local xml_file="$1"
    local entity_name="$2"
    grep -oP "(?<=<!ENTITY ${entity_name} \").*?(?=\">)" "$xml_file"
}

check_specific_job() {
    local xml_config_file="$1"
    local stage="$2"

    echo "-- Checking job status for XML: ${xml_config_file} and Stage: ${stage}..."

    local log_dir_template=$(grep -oP '(?<=<logdir>).*?(?=</logdir>)' "$xml_config_file")
    local release=$(get_xml_entity "$xml_config_file" "release")
    local project_name=$(get_xml_entity "$xml_config_file" "name")

    local log_dir=$(echo "$log_dir_template" | sed "s/&release;/${release}/" | sed "s/&name;/${project_name}/")
    local log_file="${log_dir}/jobids.list"

    if [[ ! -f "$log_file" ]]; then
        echo "Error: Log file not found: ${log_file}."
        return 1
    fi

    echo "-- Reading job ID from log file: ${log_file}..."
    local job_id=$(head -n 1 "$log_file")

    if [[ -n "$job_id" ]]; then
        echo "-- Found job ID: ${job_id}"
        local campaign_id=$(echo "$job_id" | cut -d. -f1)
        echo "-- Campaign Link: https://landscape.fnal.gov/lens/view/campaign/${campaign_id}"
        echo "-- Job Link: https://landscape.fnal.gov/lens/view/job/${job_id}"
    else
        echo "Error: No job ID found in ${log_file}."
        return 1
    fi

    echo "-- Specific job check complete!"
    return 0
}

check_all_jobs() {
    echo "-- Checking all jobs for user '${USER}'..."
    local job_queue=$(jobsub_q "$USER")

    echo "-- Job queue output:"
    echo "${job_queue}"

    echo "-- Generating campaign monitoring links..."
    echo "${job_queue}" | awk '/@/ {print $1}' | cut -d. -f1 | sort -u | while read -r campaign; do
        echo "-- Link: https://landscape.fnal.gov/lens/view/campaign/${campaign}"
    done

    echo "-- General job check complete!"
    return 0
}

USER=$(echo $USER)

if [[ "$#" -eq 2 ]]; then
    check_specific_job "$1" "$2"
elif [[ "$#" -eq 0 ]]; then
    check_all_jobs
else
    echo "Usage: source check_jobs.sh [xml_config_file] [stage]"
    echo "   or: source check_jobs.sh"
    exit 1
fi

echo "-- Job status check complete!"
