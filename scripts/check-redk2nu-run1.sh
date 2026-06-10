#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  check-redk2nu-run1.sh [--subset <all|core|fhc|rhc|detvars>] [--files <n>]
                        [--match <substring>] [--output-base <dir>] [--skip-setup]

Runs the standalone ReDk2Nu dev wrappers over the known Run 1 MC SAM surfaces
and records whether each sample produces a successful local output.

Subsets:
  all      Run everything below. Default.
  core     Run Run 1 FHC beam/dirt/current beam plus Run 1 RHC nominal samples.
  fhc      Run Run 1 FHC core samples plus generic FHC detector variations.
  rhc      Run only the Run 1 RHC samples.
  detvars  Run only the generic Run 1 FHC detector variations.

Options:
  --subset <name>       Sample subset to test. Default: all
  --files <n>           Number of files to resolve from each SAM def. Default: 1
  --match <substring>   Only run samples whose horn/group/label contains this
                        substring, e.g. --match detvar/ly or --match wiremodx
  --output-base <dir>   Base directory for logs and .local.sh outputs
  --skip-setup          Do not source .setup.sh first
  -h, --help            Show this message

Notes:
  - This script is intended to run on a uboone GPVM or equivalent environment
    where CVMFS, samweb, and lar are available.
  - EXT and data samples are intentionally excluded because ReDk2Nu only
    applies to MC inputs with dk2nu ancestry.
USAGE
}

subset="all"
files=1
match=""
skip_setup=false
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
output_base="${repo_root}/out/redk2nu_smoke"
local_runner="${repo_root}/.local.sh"

require_value() {
  if [[ $# -lt 2 ]]; then
    echo "Error: $1 requires a value." >&2
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --subset)
      require_value "$@"
      subset="$2"
      shift 2
      ;;
    --files)
      require_value "$@"
      files="$2"
      shift 2
      ;;
    --match)
      require_value "$@"
      match="$2"
      shift 2
      ;;
    --output-base)
      require_value "$@"
      output_base="$2"
      shift 2
      ;;
    --skip-setup)
      skip_setup=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'." >&2
      usage >&2
      exit 1
      ;;
  esac
done

case "${subset}" in
  all|core|fhc|rhc|detvars)
    ;;
  *)
    echo "Error: unsupported subset '${subset}'." >&2
    exit 1
    ;;
esac

if ! [[ "${files}" =~ ^[0-9]+$ ]] || (( files <= 0 )); then
  echo "Error: --files must be a positive integer." >&2
  exit 1
fi

if [[ ! -x "${local_runner}" ]]; then
  echo "Error: local runner not found: ${local_runner}" >&2
  exit 1
fi

if [[ "${skip_setup}" != true ]]; then
  # UPS/CVMFS setup scripts are not strict-mode clean; relax while sourcing.
  set +eu
  # shellcheck disable=SC1091
  source "${repo_root}/.setup.sh"
  set -eu
fi

if ! command -v samweb >/dev/null 2>&1; then
  echo "Error: samweb is not available. Source the uboone products setup first." >&2
  exit 1
fi

if ! command -v lar >/dev/null 2>&1; then
  echo "Error: lar is not available. Source the uboone products setup first." >&2
  exit 1
fi

slugify() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | tr -cs 'a-z0-9' '_' \
    | sed 's/^_*//; s/_*$//'
}

reason_from_log() {
  local log_file="$1"
  local reason=""

  if grep -Fq 'Found 0 allowed dk2nu entries' "${log_file}"; then
    echo 'Found 0 allowed dk2nu entries'
    return
  fi

  reason="$(sed -n 's/^Error: \(.*\)$/\1/p' "${log_file}" | tail -n 1)"
  if [[ -n "${reason}" ]]; then
    printf '%s\n' "${reason}"
    return
  fi

  reason="$(sed -n 's/^Art has completed and will exit with status \(.*\)\.$/art exit status \1/p' "${log_file}" | tail -n 1)"
  if [[ -n "${reason}" ]]; then
    printf '%s\n' "${reason}"
    return
  fi

  tail -n 1 "${log_file}" | tr '\t' ' ' | tr -s ' '
}

declare -a samples=()

add_sample() {
  local horn="$1"
  local sample_group="$2"
  local label="$3"
  local samdef="$4"
  local fhicl_file="dev/run_stage_redk2nu_dev.fcl"
  if [[ "${horn}" == "rhc" ]]; then
    fhicl_file="dev/run_stage_redk2nu_rhc_dev.fcl"
  fi
  samples+=("${horn}|${sample_group}|${label}|${fhicl_file}|${samdef}")
}

add_fhc_core() {
  add_sample "fhc" "beam" "legacy_sample0" "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0"
  add_sample "fhc" "beam" "legacy_sample1" "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1"
  add_sample "fhc" "beam" "legacy_sample2" "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2"
  add_sample "fhc" "beam" "legacy_sample3" "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3"
  add_sample "fhc" "beam" "merged_pandora" "New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2"
  add_sample "fhc" "dirt" "legacy_dirt_sample0" "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0"
  add_sample "fhc" "dirt" "legacy_dirt_sample1" "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1"
  add_sample "fhc" "strangeness" "dedicated_strangeness" "prod_strange_resample_fhc_run1_fhc_reco2_reco2"
}

add_fhc_detvars() {
  add_sample "fhc" "detvar" "cv" "prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"
  add_sample "fhc" "detvar" "ly_attenuation" "prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2"
  add_sample "fhc" "detvar" "ly_rayleigh" "prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2"
  add_sample "fhc" "detvar" "ly_down" "prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2"
  add_sample "fhc" "detvar" "sce" "prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2"
  add_sample "fhc" "detvar" "recomb2" "prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2"
  add_sample "fhc" "detvar" "wiremodx" "prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2"
  add_sample "fhc" "detvar" "wiremodyz" "prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2"
  add_sample "fhc" "detvar" "wiremodthetaxz" "prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2"
  add_sample "fhc" "detvar" "wiremodthetayz" "prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2"
}

add_rhc_core() {
  add_sample "rhc" "beam" "nu_overlay" "prodgenie_numi_rhc_nu_overlay_v08_00_00_54_run1_reco2_reco2"
  add_sample "rhc" "beam" "intrinsic_nue" "prodgenie_numi_rhc_intrinsic_nue_overlay_v08_00_00_54_run1_reco2_reco2_reco2"
  add_sample "rhc" "dirt" "dirt_overlay" "prod_extunbiased_numi_rhc_dirt_overlay_run1_reco2_v08_00_00_67_reco2"
}

case "${subset}" in
  all)
    add_fhc_core
    add_fhc_detvars
    add_rhc_core
    ;;
  core)
    add_fhc_core
    add_rhc_core
    ;;
  fhc)
    add_fhc_core
    add_fhc_detvars
    ;;
  rhc)
    add_rhc_core
    ;;
  detvars)
    add_fhc_detvars
    ;;
esac

if [[ -n "${match}" ]]; then
  declare -a matched=()
  for entry in "${samples[@]}"; do
    IFS='|' read -r horn sample_group label _ _ <<< "${entry}"
    if [[ "${horn}/${sample_group}/${label}" == *"${match}"* ]]; then
      matched+=("${entry}")
    fi
  done
  if (( ${#matched[@]} == 0 )); then
    echo "Error: no samples match '--match ${match}'." >&2
    exit 1
  fi
  samples=("${matched[@]}")
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
run_dir="${output_base}/run1_redk2nu_${subset}_${timestamp}"
logs_dir="${run_dir}/logs"
mkdir -p "${logs_dir}"

summary_tsv="${run_dir}/summary.tsv"
printf 'index\thorn\tgroup\tlabel\tstatus\tsamdef\tfhicl\treason\tlog\n' > "${summary_tsv}"

trap 'echo; echo "Interrupted. Partial summary: ${summary_tsv}" >&2; exit 130' INT TERM

pass_count=0
fail_count=0
total_count="${#samples[@]}"

echo "Writing logs and outputs under: ${run_dir}"
echo "Testing ${total_count} Run 1 sample(s) with ${files} file(s) each."

for idx in "${!samples[@]}"; do
  IFS='|' read -r horn sample_group label fhicl_file samdef <<< "${samples[$idx]}"
  sample_slug="$(slugify "${horn}_${sample_group}_${label}")"
  log_file="${logs_dir}/$(printf '%02d' "$((idx + 1))")_${sample_slug}.log"
  status="FAIL"
  reason=""

  echo
  echo "[$((idx + 1))/${total_count}] ${horn}/${sample_group}/${label}"
  echo "  SAM_DEF=${samdef}"
  echo "  FHiCL=${fhicl_file}"

  cmd_status=0
  SAM_DEF="${samdef}" OUTPUT_BASE_DIR="${run_dir}/out" \
    bash "${local_runner}" "${fhicl_file}" "${files}" >"${log_file}" 2>&1 || cmd_status=$?

  if [[ ${cmd_status} -eq 0 ]]; then
    status="PASS"
    reason="ok"
    pass_count=$((pass_count + 1))
  else
    reason="$(reason_from_log "${log_file}")"
    fail_count=$((fail_count + 1))
  fi

  reason="${reason//$'\t'/ }"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$((idx + 1))" \
    "${horn}" \
    "${sample_group}" \
    "${label}" \
    "${status}" \
    "${samdef}" \
    "${fhicl_file}" \
    "${reason}" \
    "${log_file}" >> "${summary_tsv}"

  echo "  ${status}: ${reason}"
done

echo
cut -f1-5,8 "${summary_tsv}" | column -t -s$'\t' || cat "${summary_tsv}"
echo
echo "Summary written to: ${summary_tsv}"
echo "Passed: ${pass_count}"
echo "Failed: ${fail_count}"

if (( fail_count > 0 )); then
  exit 1
fi
