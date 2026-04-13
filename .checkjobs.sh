#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  ./.checkjobs.sh
  ./.checkjobs.sh <xml_config_file> <stage>

With no arguments, print the current jobsub queue for $USER and Lens links.
With an XML and stage, print the job IDs recorded for that campaign stage.
USAGE
}

get_xml_entity() {
    local xml_file="$1"
    local entity_name="$2"

    awk -v entity="${entity_name}" '
        $0 ~ "<!ENTITY[[:space:]]+" entity "[[:space:]]+\"" {
            line = $0
            sub(".*<!ENTITY[[:space:]]+" entity "[[:space:]]+\"", "", line)
            sub("\">.*", "", line)
            print line
            exit
        }
    ' "${xml_file}"
}

get_stage_logdir() {
    local xml_file="$1"
    local stage="$2"

    awk -v target="${stage}" '
        $0 ~ "<stage[[:space:]]+name=\"" target "\"" {
            in_stage = 1
        }
        in_stage && /<logdir>/ {
            line = $0
            sub(".*<logdir>", "", line)
            sub("</logdir>.*", "", line)
            print line
            exit
        }
        in_stage && /<\/stage>/ {
            in_stage = 0
        }
    ' "${xml_file}"
}

expand_known_entities() {
    local xml_file="$1"
    local value="$2"
    local release user project_name

    release="$(get_xml_entity "${xml_file}" "release")"
    user="$(get_xml_entity "${xml_file}" "user")"
    project_name="$(get_xml_entity "${xml_file}" "name")"

    value="${value//&release;/${release}}"
    value="${value//&user;/${user}}"
    value="${value//&name;/${project_name}}"

    printf '%s\n' "${value}"
}

print_lens_links() {
    awk 'NF {print "-- Link: https://landscape.fnal.gov/lens/view/job/" $1}'
}

check_specific_job() {
    local xml_config_file="$1"
    local stage="$2"
    local log_dir_template log_dir log_file

    if [[ ! -f "${xml_config_file}" ]]; then
        echo "Error: XML file not found: ${xml_config_file}" >&2
        return 1
    fi

    log_dir_template="$(get_stage_logdir "${xml_config_file}" "${stage}")"
    if [[ -z "${log_dir_template}" ]]; then
        echo "Error: stage '${stage}' not found in ${xml_config_file}, or it has no <logdir>." >&2
        return 1
    fi

    log_dir="$(expand_known_entities "${xml_config_file}" "${log_dir_template}")"
    log_file="${log_dir}/jobids.list"

    echo "-- XML: ${xml_config_file}"
    echo "-- Stage: ${stage}"
    echo "-- Log dir: ${log_dir}"

    if [[ ! -f "${log_file}" ]]; then
        echo "Error: job ID file not found: ${log_file}" >&2
        echo "Run project.py --xml ${xml_config_file} --stage ${stage} --check, or confirm the stage was submitted." >&2
        return 1
    fi

    echo "-- Job IDs:"
    cat "${log_file}"
    echo "-- Lens links:"
    print_lens_links < "${log_file}"
}

check_all_jobs() {
    if ! command -v jobsub_q >/dev/null 2>&1; then
        echo "Error: jobsub_q is not available in PATH." >&2
        return 1
    fi

    echo "-- Checking all jobs for user '${USER}'..."
    jobsub_q --user="${USER}" | tee /tmp/checkjobs.$$.queue
    echo "-- Lens links:"
    awk '/@/ {print $1}' /tmp/checkjobs.$$.queue | print_lens_links
    rm -f /tmp/checkjobs.$$.queue
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
elif [[ "$#" -eq 2 ]]; then
    check_specific_job "$1" "$2"
elif [[ "$#" -eq 0 ]]; then
    check_all_jobs
else
    usage >&2
    exit 1
fi
