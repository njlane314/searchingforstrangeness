#!/usr/bin/env bash

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "${REPO_DIR}"

source "runtime/scripts/initsrc.sh"

usage() {
    cat <<EOF
Usage:
  $(basename "$0") <fhicl_file> <num_files>
  $(basename "$0") <fhicl_file> <input.root>

Environment:
  SAM_DEF         SAM definition used when <num_files> is provided.
  OUTPUT_BASE_DIR Directory for local outputs.
EOF
}

is_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

abs_path() {
    local target="$1"
    local target_dir
    target_dir="$(cd "$(dirname "$target")" && pwd -P)"
    printf '%s/%s\n' "${target_dir}" "$(basename "$target")"
}

fhicl_has_tfile_output() {
    grep -Eq 'TFileService|services\.TFileService\.fileName' "${FHICL_FILE}"
}

fhicl_has_art_output() {
    grep -Eq 'outputs\.out1\.fileName|module_type:[[:space:]]*"?RootOutput"?' "${FHICL_FILE}"
}

fetch_files_from_sam() {
    mapfile -t FETCHED_FILES < <(samweb list-files defname:"${SAM_DEF}" | head -n "${NUM_FILES}")

    if [ "${#FETCHED_FILES[@]}" -eq 0 ]; then
        echo "Error: No files found in SAM definition '${SAM_DEF}'." >&2
        exit 1
    fi

    echo "Fetched ${#FETCHED_FILES[@]} file(s) from '${SAM_DEF}'." >&2
}

resolve_sam_file_path() {
    local file_name="$1"
    local file_dir

    file_dir="$(samweb locate-file "${file_name}" | grep -o '/pnfs/.*' | sed 's/(\([0-9]*@[a-z0-9]*\))//g' | head -n 1)"
    if [ -z "${file_dir}" ]; then
        echo "Error: Could not locate directory for file '${file_name}'." >&2
        return 1
    fi

    printf '%s/%s\n' "${file_dir}" "${file_name}"
}

run_lar() {
    local input_path="$1"
    local stem="$2"
    local lar_succeeded=0

    local -a cmd
    cmd=(lar -c "${FHICL_FILE}" -s "${input_path}")

    if [ "${HAS_ART_OUTPUT}" -eq 1 ]; then
        CURRENT_EVENT_OUTPUT="${RUN_DIR}/${stem}_events.root"
        cmd+=(-o "${CURRENT_EVENT_OUTPUT}")
    else
        CURRENT_EVENT_OUTPUT=""
    fi

    if [ "${HAS_TFILE_OUTPUT}" -eq 1 ]; then
        CURRENT_HIST_OUTPUT="${RUN_DIR}/${stem}_hist.root"
        cmd+=(-T "${CURRENT_HIST_OUTPUT}")
    else
        CURRENT_HIST_OUTPUT=""
    fi

    printf 'Running:'
    printf ' %q' "${cmd[@]}"
    printf '\n'

    if "${cmd[@]}"; then
        lar_succeeded=1
    fi

    if [ "${lar_succeeded}" -ne 1 ]; then
        echo "Error: FHiCL processing failed for '${input_path}'." >&2
        rm -f "${CURRENT_EVENT_OUTPUT:-}" "${CURRENT_HIST_OUTPUT:-}"
        return 1
    fi

    if [ -n "${CURRENT_EVENT_OUTPUT}" ] && [ -f "${CURRENT_EVENT_OUTPUT}" ]; then
        SUCCESSFUL_EVENT_OUTPUTS+=("${CURRENT_EVENT_OUTPUT}")
    elif [ -n "${CURRENT_EVENT_OUTPUT}" ]; then
        echo "Warning: Expected art output '${CURRENT_EVENT_OUTPUT}' was not created." >&2
    fi

    if [ -n "${CURRENT_HIST_OUTPUT}" ] && [ -f "${CURRENT_HIST_OUTPUT}" ]; then
        SUCCESSFUL_HIST_OUTPUTS+=("${CURRENT_HIST_OUTPUT}")
    elif [ -n "${CURRENT_HIST_OUTPUT}" ]; then
        echo "Warning: Expected histogram output '${CURRENT_HIST_OUTPUT}' was not created." >&2
    fi
}

process_sam_inputs() {
    local counter=0
    local file_name
    local input_path
    local stem

    for file_name in "${FETCHED_FILES[@]}"; do
        echo "Processing SAM file: ${file_name}"
        input_path="$(resolve_sam_file_path "${file_name}")" || continue
        if [ ! -f "${input_path}" ]; then
            echo "Error: File not found at '${input_path}'. Skipping..." >&2
            continue
        fi

        stem="$(printf '%03d_%s' "${counter}" "$(basename "${file_name}" .root)")"
        run_lar "${input_path}" "${stem}" || true
        counter=$((counter + 1))
    done
}

process_local_input() {
    local input_path="$1"
    local stem

    if [ ! -f "${input_path}" ]; then
        echo "Error: Local input '${input_path}' does not exist." >&2
        exit 1
    fi

    input_path="$(abs_path "${input_path}")"
    stem="$(basename "${input_path}" .root)"
    echo "Processing local file: ${input_path}"
    run_lar "${input_path}" "${stem}"
}

finalize_hist_outputs() {
    if [ "${#SUCCESSFUL_HIST_OUTPUTS[@]}" -eq 0 ]; then
        return
    fi

    COMBINED_HIST_OUTPUT="${RUN_DIR}/${RUN_STEM}_hist.root"
    if [ "${#SUCCESSFUL_HIST_OUTPUTS[@]}" -eq 1 ]; then
        cp -f "${SUCCESSFUL_HIST_OUTPUTS[0]}" "${COMBINED_HIST_OUTPUT}"
    else
        hadd -f "${COMBINED_HIST_OUTPUT}" "${SUCCESSFUL_HIST_OUTPUTS[@]}"
    fi

    echo "Combined histogram output: ${COMBINED_HIST_OUTPUT}"
}

finalize_event_outputs() {
    if [ "${#SUCCESSFUL_EVENT_OUTPUTS[@]}" -eq 0 ]; then
        return
    fi

    EVENT_OUTPUT_LIST="${RUN_DIR}/${RUN_STEM}_events.list"
    printf '%s\n' "${SUCCESSFUL_EVENT_OUTPUTS[@]}" > "${EVENT_OUTPUT_LIST}"

    echo "Event-output list: ${EVENT_OUTPUT_LIST}"
    if [ "${#SUCCESSFUL_EVENT_OUTPUTS[@]}" -eq 1 ]; then
        echo "Single event output: ${SUCCESSFUL_EVENT_OUTPUTS[0]}"
    fi
}

if [ "$#" -ne 2 ]; then
    usage
    exit 1
fi

FHICL_FILE="$1"
INPUT_SPEC="$2"

if [ ! -f "${FHICL_FILE}" ]; then
    echo "Error: FHiCL file '${FHICL_FILE}' does not exist." >&2
    exit 1
fi

FHICL_FILE="$(abs_path "${FHICL_FILE}")"
FHICL_BASE="$(basename "${FHICL_FILE}" .fcl | sed 's/^run_//')"
SAM_DEF="${SAM_DEF:-prod_strange_resample_fhc_run2_fhc_reco2_reco2}"
OUTPUT_BASE_DIR="${OUTPUT_BASE_DIR:-${REPO_DIR}/local-output}"
RUN_STEM="${SAM_DEF}_${FHICL_BASE}"

if is_integer "${INPUT_SPEC}"; then
    NUM_FILES="${INPUT_SPEC}"
    RUN_STEM="${RUN_STEM}_${NUM_FILES}"
else
    NUM_FILES=1
    RUN_STEM="$(basename "${INPUT_SPEC}" .root)_${FHICL_BASE}_local"
fi

RUN_DIR="${OUTPUT_BASE_DIR}/${RUN_STEM}"
if [ -e "${RUN_DIR}" ]; then
    RUN_DIR="${RUN_DIR}_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "${RUN_DIR}"

SUCCESSFUL_EVENT_OUTPUTS=()
SUCCESSFUL_HIST_OUTPUTS=()
HAS_ART_OUTPUT=0
HAS_TFILE_OUTPUT=0

if fhicl_has_art_output; then
    HAS_ART_OUTPUT=1
fi

if fhicl_has_tfile_output; then
    HAS_TFILE_OUTPUT=1
fi

echo "Writing outputs under: ${RUN_DIR}"

if is_integer "${INPUT_SPEC}"; then
    echo "Starting SAM-backed local test for '${SAM_DEF}' with '${FHICL_FILE}'."
    fetch_files_from_sam
    process_sam_inputs
else
    echo "Starting local-file test with '${FHICL_FILE}'."
    process_local_input "${INPUT_SPEC}"
fi

if [ "${#SUCCESSFUL_EVENT_OUTPUTS[@]}" -eq 0 ] && [ "${#SUCCESSFUL_HIST_OUTPUTS[@]}" -eq 0 ]; then
    echo "Error: No successful outputs were produced." >&2
    exit 1
fi

finalize_hist_outputs
finalize_event_outputs

echo "Process complete!"
