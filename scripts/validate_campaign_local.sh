#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  validate_campaign_local.sh [--workflow <mc|data|amarantin|fullchain>] [--samdef <def>] [--files <n>]
  validate_campaign_local.sh [--workflow <mc|data|amarantin|fullchain>] --input <input.root>

Runs the checked-in dev FHiCL wrappers locally through .local.sh so you can
validate the campaign path before submitting to the grid.

Workflows:
  mc         Run staged MC validation: evtw -> image -> sel
             Default evtw config is cv, matching the active campaign XMLs.
  data       Run staged data/EXT-style validation: image -> sel_data
  amarantin  Run the compact downstream ntuple surface for amarantin
  fullchain  Run the single-process dev fullchain wrapper

Options:
  --workflow <name>       Validation workflow. Default: mc
  --samdef <def>          Source SAM definition when no --input is given
  --files <n>             Number of SAM files to resolve. Default: 1
  --input <input.root>    Local ROOT file to use instead of SAM
  --evtw-config <name>    One of: cv, extragenie1, extragenie2, extragenie3,
                          extragenie4, extragenie5. Default: cv
  --output-base <dir>     Base directory for local-output runs

Notes:
  - staged mc/data validation requires exactly one starting input file when
    using SAM, because local stage handoff is one ROOT file at a time
  - fullchain is convenient for debugging, but it is not the same surface as
    the staged campaign path
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
local_runner="${repo_root}/.local.sh"

workflow="mc"
samdef=""
files=1
input_path=""
evtw_config="cv"
output_base="${repo_root}/local-output/validation"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workflow)
      if [[ $# -lt 2 ]]; then
        echo "Error: --workflow expects a value." >&2
        exit 1
      fi
      workflow="$2"
      shift 2
      ;;
    --samdef)
      if [[ $# -lt 2 ]]; then
        echo "Error: --samdef expects a value." >&2
        exit 1
      fi
      samdef="$2"
      shift 2
      ;;
    --files)
      if [[ $# -lt 2 ]]; then
        echo "Error: --files expects an integer." >&2
        exit 1
      fi
      files="$2"
      shift 2
      ;;
    --input)
      if [[ $# -lt 2 ]]; then
        echo "Error: --input expects a path." >&2
        exit 1
      fi
      input_path="$2"
      shift 2
      ;;
    --evtw-config)
      if [[ $# -lt 2 ]]; then
        echo "Error: --evtw-config expects a value." >&2
        exit 1
      fi
      evtw_config="$2"
      shift 2
      ;;
    --output-base)
      if [[ $# -lt 2 ]]; then
        echo "Error: --output-base expects a path." >&2
        exit 1
      fi
      output_base="$2"
      shift 2
      ;;
    *)
      echo "Error: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! -f "${local_runner}" ]]; then
  echo "Error: local runner not found: ${local_runner}" >&2
  exit 1
fi

if [[ -n "${input_path}" && -n "${samdef}" ]]; then
  echo "Error: use either --input or --samdef, not both." >&2
  exit 1
fi

if ! [[ ${files} =~ ^[0-9]+$ ]] || (( files <= 0 )); then
  echo "Error: --files must be a positive integer." >&2
  exit 1
fi

case "${workflow}" in
  mc|data|amarantin|fullchain)
    ;;
  *)
    echo "Error: unsupported workflow: ${workflow}" >&2
    exit 1
    ;;
esac

case "${evtw_config}" in
  cv|extragenie1|extragenie2|extragenie3|extragenie4|extragenie5)
    ;;
  *)
    echo "Error: unsupported --evtw-config value: ${evtw_config}" >&2
    exit 1
    ;;
esac

default_mc_samdef="prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0"
default_data_samdef="prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2"

if [[ -z "${input_path}" && -z "${samdef}" ]]; then
  case "${workflow}" in
    mc|amarantin|fullchain)
      samdef="${default_mc_samdef}"
      ;;
    data)
      samdef="${default_data_samdef}"
      ;;
  esac
fi

if [[ -z "${input_path}" && ( "${workflow}" == "mc" || "${workflow}" == "data" ) && ${files} -ne 1 ]]; then
  echo "Error: staged ${workflow} validation requires --files 1 so the local handoff stays on a single ROOT file." >&2
  exit 1
fi

if [[ -n "${input_path}" && ! -f "${input_path}" ]]; then
  echo "Error: input file not found: ${input_path}" >&2
  exit 1
fi

run_step() {
  local fhicl_file="$1"
  local input_spec="$2"
  local log_file

  log_file="$(mktemp)"
  LAST_EVENT_OUTPUT=""
  LAST_EVENT_LIST=""
  LAST_HIST_OUTPUT=""
  LAST_RUN_DIR=""
  trap 'rm -f "${log_file}"' RETURN

  echo "Running local validation step: ${fhicl_file}"

  if [[ -n "${samdef}" && "${input_spec}" =~ ^[0-9]+$ ]]; then
    SAM_DEF="${samdef}" OUTPUT_BASE_DIR="${output_base}" bash "${local_runner}" "${fhicl_file}" "${input_spec}" | tee "${log_file}"
  else
    OUTPUT_BASE_DIR="${output_base}" bash "${local_runner}" "${fhicl_file}" "${input_spec}" | tee "${log_file}"
  fi

  LAST_EVENT_OUTPUT="$(sed -n 's/^Single event output: //p' "${log_file}" | tail -n 1)"
  LAST_EVENT_LIST="$(sed -n 's/^Event-output list: //p' "${log_file}" | tail -n 1)"
  LAST_HIST_OUTPUT="$(sed -n 's/^Combined histogram output: //p' "${log_file}" | tail -n 1)"
  LAST_RUN_DIR="$(sed -n 's/^Writing outputs under: //p' "${log_file}" | tail -n 1)"

  if [[ -z "${LAST_EVENT_OUTPUT}" && -n "${LAST_EVENT_LIST}" && -f "${LAST_EVENT_LIST}" ]]; then
    if [[ $(wc -l < "${LAST_EVENT_LIST}") -eq 1 ]]; then
      LAST_EVENT_OUTPUT="$(head -n 1 "${LAST_EVENT_LIST}")"
    fi
  fi
}

require_single_event_output() {
  local stage_label="$1"

  if [[ -z "${LAST_EVENT_OUTPUT}" ]]; then
    echo "Error: ${stage_label} did not produce a single event ROOT file for the next stage." >&2
    if [[ -n "${LAST_EVENT_LIST}" ]]; then
      echo "Event list: ${LAST_EVENT_LIST}" >&2
    fi
    exit 1
  fi
}

print_summary() {
  echo
  echo "Local validation summary:"
  if [[ -n "${LAST_RUN_DIR}" ]]; then
    echo "  final run dir : ${LAST_RUN_DIR}"
  fi
  if [[ -n "${LAST_HIST_OUTPUT}" ]]; then
    echo "  hist output   : ${LAST_HIST_OUTPUT}"
  fi
  if [[ -n "${LAST_EVENT_OUTPUT}" ]]; then
    echo "  event output  : ${LAST_EVENT_OUTPUT}"
  elif [[ -n "${LAST_EVENT_LIST}" ]]; then
    echo "  event list    : ${LAST_EVENT_LIST}"
  fi
}

case "${workflow}" in
  mc)
    case "${evtw_config}" in
      cv) evtw_fhicl="dev/run_stage_evtw_cv_dev.fcl" ;;
      extragenie1) evtw_fhicl="dev/run_stage_evtw_extragenie1_dev.fcl" ;;
      extragenie2) evtw_fhicl="dev/run_stage_evtw_extragenie2_dev.fcl" ;;
      extragenie3) evtw_fhicl="dev/run_stage_evtw_extragenie3_dev.fcl" ;;
      extragenie4) evtw_fhicl="dev/run_stage_evtw_extragenie4_dev.fcl" ;;
      extragenie5) evtw_fhicl="dev/run_stage_evtw_extragenie5_dev.fcl" ;;
    esac

    if [[ -n "${input_path}" ]]; then
      first_input="${input_path}"
    else
      first_input="${files}"
    fi

    run_step "${evtw_fhicl}" "${first_input}"
    require_single_event_output "eventweight stage"
    run_step "dev/run_stage_image_dev.fcl" "${LAST_EVENT_OUTPUT}"
    require_single_event_output "image stage"
    run_step "dev/run_stage_sel_dev.fcl" "${LAST_EVENT_OUTPUT}"
    print_summary
    ;;
  data)
    if [[ -n "${input_path}" ]]; then
      first_input="${input_path}"
    else
      first_input="${files}"
    fi

    run_step "dev/run_stage_image_dev.fcl" "${first_input}"
    require_single_event_output "image stage"
    run_step "dev/run_stage_sel_data_dev.fcl" "${LAST_EVENT_OUTPUT}"
    print_summary
    ;;
  amarantin)
    if [[ -n "${input_path}" ]]; then
      first_input="${input_path}"
    else
      first_input="${files}"
    fi

    run_step "dev/amarantin_local_dev.fcl" "${first_input}"
    print_summary
    ;;
  fullchain)
    if [[ -n "${input_path}" ]]; then
      first_input="${input_path}"
    else
      first_input="${files}"
    fi

    run_step "dev/run_stage_fullchain_dev.fcl" "${first_input}"
    print_summary
    ;;
esac
