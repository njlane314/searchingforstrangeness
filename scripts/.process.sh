#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

cleanup_temp_dir() {
    if [[ -n "${TEMP_DIR:-}" && -d "${TEMP_DIR}" ]]; then
        echo "Cleaning up temporary files in ${TEMP_DIR}..."
        rm -rf "${TEMP_DIR}" || {
            echo "Warning: Failed to remove temporary directory '${TEMP_DIR}'. Manual cleanup might be needed."
        }
        echo "Temporary files cleaned."
    fi
}

trap cleanup_temp_dir EXIT

check_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "Error: required command '$1' not found. Exiting..."
        exit 1
    }
}

fetch_files_from_sam() {
    local files_list=$(samweb list-files defname:"${SAM_DEF}" | head -n "${NUM_FILES}")

    if [ -z "${files_list}" ]; then
        echo "Error: No files found in SAM definition '${SAM_DEF}'. Exiting..."
        exit 1
    fi
    echo "Files fetched successfully."
    echo "${files_list}"
}

process_files_with_lar() {
    local files="$1"
    mkdir -p "${TEMP_DIR}"

    local counter=0
    for file in ${files}; do
        echo "Processing file: ${file}..."

        local filedir=$(samweb locate-file "${file}" | grep -o '/pnfs/.*' | sed 's/(\([0-9]*@[a-z0-9]*\))//g' | head -n 1)
        if [ -z "${filedir}" ]; then
            echo "Error: Could not locate directory for file: ${file}. Skipping..."
            continue
        fi

        local filepath="${filedir}/${file}"
        if [ ! -f "${filepath}" ]; then
            echo "Error: File not found at: ${filepath}. Skipping..."
            continue
        fi

        local outputfile="${TEMP_DIR}/output_${counter}.root"
        echo "Running: lar -c ${FHICL_FILE} -s ${filepath} -T ${outputfile}"
        lar -c "${FHICL_FILE}" -s "${filepath}" -T "${outputfile}" || {
            echo "Error: FHiCL processing failed for file: ${file}."
            continue
        }
        echo "FHiCL processing successful. Output: ${outputfile}"
        counter=$((counter + 1))
    done
}

combine_output_files() {
    echo "Combining all the output ROOT files..."
    local outputfiles=$(find "${TEMP_DIR}" -maxdepth 1 -name "*.root")

    if [ -z "${outputfiles}" ]; then
        echo "Error: No output ROOT files found to combine! Exiting..."
        exit 1
    fi

    echo "Combining files into: ${COMBINED_OUTPUT}"
    hadd -f "${COMBINED_OUTPUT}" ${outputfiles} || {
        echo "Error: Combining ROOT files failed!"
        exit 1
    }
    echo "Successfully combined ROOT files into: ${COMBINED_OUTPUT}"
}

FHICL_FILE="$1"
NUM_FILES="$2"
SAM_DEF="New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2"
OUTPUT_BASE_DIR="/exp/uboone/data/users/$USER/analysis"
FHICL_BASE=$(basename "${FHICL_FILE}" .fcl | sed 's/^run_//')
COMBINED_OUTPUT="${OUTPUT_BASE_DIR}/${SAM_DEF}_${FHICL_BASE}_${NUM_FILES}_new_analysis.root"
TEMP_DIR="${OUTPUT_BASE_DIR}/temp_root_files"

check_command samweb
check_command lar
check_command hadd

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <fhicl_file> <num_files>"
    exit 1
fi

echo "Starting process for SAM definition: ${SAM_DEF}"

FETCHED_FILES=$(fetch_files_from_sam)

process_files_with_lar "${FETCHED_FILES}"

combine_output_files

echo "Process complete!"

