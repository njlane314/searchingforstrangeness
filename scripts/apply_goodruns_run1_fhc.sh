#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  apply_goodruns_run1_fhc.sh [--condition <expr>] [--dry-run]

Applies the good-runs condition to the Run 1 NuMI FHC beam, generic detvar,
dirt, and EXT source definitions listed in scripts/README.md, generating
consistent output names. Dedicated strangeness Reco2 SAM definitions are
excluded because they were already produced from good-run-filtered inputs.

Options:
  --condition <expr>  Override the condition (default: "goodruns: 1").
  --dry-run           Print source -> output mappings without creating defs.
  -h, --help          Show this help message.
USAGE
}

condition="goodruns: 1"
dry_run=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --condition)
      condition="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'" >&2
      usage >&2
      exit 1
      ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
apply_script="${script_dir}/apply_goodruns.sh"

if [[ ! -x "${apply_script}" ]]; then
  echo "Error: ${apply_script} is not executable." >&2
  exit 1
fi

sources=(
  # Beam
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3"

  # Detector variations (Run 1 FHC)
  "prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2"
  "prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2"

  # Dirt
  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0"
  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1"

  # EXT & Data (Run 1 FHC)
  "prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2"
)

goodruns_name() {
  local source="$1"

  declare -A short_names=(
    ["prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0"]="nl_run1_fhc_beam_sample0_goodruns"
    ["prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1"]="nl_run1_fhc_beam_sample1_goodruns"
    ["prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2"]="nl_run1_fhc_beam_sample2_goodruns"
    ["prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3"]="nl_run1_fhc_beam_sample3_goodruns"
    ["prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"]="nl_run1_fhc_detvar_cv_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_ly_supp75_att8m_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_ly_rayleigh_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_ly_down_goodruns"
    ["prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2"]="nl_run1_fhc_detvar_sce_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_recomb2_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2"]="nl_run1_fhc_detvar_wiremodx_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_wiremodyz_goodruns"
    ["prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2"]="nl_run1_fhc_detvar_wiremodthetaxz_goodruns"
    ["prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2"]="nl_run1_fhc_detvar_wiremodthetayz_goodruns"
    ["prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0"]="nl_run1_fhc_dirt_sample0_goodruns"
    ["prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1"]="nl_run1_fhc_dirt_sample1_goodruns"
    ["prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2"]="nl_run1_fhc_ext_goodruns"
  )

  if [[ -n "${short_names[${source}]:-}" ]]; then
    echo "${short_names[${source}]}"
    return
  fi

  if [[ "${source}" == *_goodruns ]]; then
    echo "Error: source '${source}' already contains a goodruns suffix." >&2
    exit 1
  fi

  if [[ "${source}" == *"_good_reco2" ]]; then
    echo "${source/_good_reco2/_goodruns_reco2}"
    return
  fi

  echo "${source}_goodruns"
}

declare -A seen_targets

for source in "${sources[@]}"; do
  target="$(goodruns_name "${source}")"

  if [[ -n "${seen_targets[${target}]:-}" ]]; then
    echo "Error: duplicate target definition '${target}'." >&2
    exit 1
  fi
  seen_targets["${target}"]=1

  if [[ "${dry_run}" == true ]]; then
    printf '%s -> %s\n' "${source}" "${target}"
    continue
  fi

  "${apply_script}" "${source}" "${target}" "${condition}"
  echo
  echo "---"
  echo
  sleep 1
done
