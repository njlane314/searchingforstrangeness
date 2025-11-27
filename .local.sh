#!/bin/bash

set -e

source "assets/initsrc.sh"

fetch_files_from_sam() {
    local files_list=$(samweb list-files defname:"${SAM_DEF}" | head -n "${NUM_FILES}")

    if [ -z "${files_list}" ]; then
        echo "Error: No files found in SAM definition '${SAM_DEF}'. Exiting..."
        exit 1
    fi
    echo "Files fetched successfully." >&2
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
        rm -f "${outputfile}"
        echo "Running: lar -c ${FHICL_FILE} -s ${filepath} -T ${outputfile}"
        if lar -c "${FHICL_FILE}" -s "${filepath}" -T "${outputfile}"; then
            echo "FHiCL processing successful. Output: ${outputfile}"
            SUCCESSFUL_OUTPUTS="${SUCCESSFUL_OUTPUTS} ${outputfile}"
            counter=$((counter + 1))
        else
            echo "Error: FHiCL processing failed for file: ${file}."
            rm -f "${outputfile}"
            continue
        fi
    done
}

combine_output_files() {
    echo "Combining all the output ROOT files..."

    if [ -z "${SUCCESSFUL_OUTPUTS}" ]; then
        echo "Error: No successful output ROOT files found to combine! Exiting..."
        exit 1
    fi

    set -- ${SUCCESSFUL_OUTPUTS}
    echo "Combining files into: ${COMBINED_OUTPUT}"
    hadd -f "${COMBINED_OUTPUT}" "$@" || {
        echo "Error: Combining ROOT files failed!"
        exit 1
    }
    echo "Successfully combined ROOT files into: ${COMBINED_OUTPUT}"
}

cleanup_temp_dir() {
    echo "Cleaning up temporary files in ${TEMP_DIR}..."
    rm -rf "${TEMP_DIR}" || {
        echo "Warning: Failed to remove temporary directory '${TEMP_DIR}'. Manual cleanup might be needed."
    }
    echo "Temporary files cleaned."
}

FHICL_FILE="$1"
NUM_FILES="$2"
#SAM_DEF="prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"
SAM_DEF="prod_strange_resample_fhc_run2_fhc_reco2_reco2"
#SAM_DEF="nl_prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2_3000"
OUTPUT_BASE_DIR="/exp/uboone/data/users/$USER/analysis"
FHICL_BASE=$(basename "${FHICL_FILE}" .fcl | sed 's/^run_//')
COMBINED_OUTPUT="${OUTPUT_BASE_DIR}/${SAM_DEF}_${FHICL_BASE}_${NUM_FILES}_new_analysis.root"
TEMP_DIR="${OUTPUT_BASE_DIR}/temp_root_files"
SUCCESSFUL_OUTPUTS=""

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <fhicl_file> <num_files>"
    exit 1
fi

echo "Starting process for SAM definition: ${SAM_DEF}"

FETCHED_FILES=$(fetch_files_from_sam)

process_files_with_lar "${FETCHED_FILES}"

combine_output_files

cleanup_temp_dir

echo "Process complete!"
