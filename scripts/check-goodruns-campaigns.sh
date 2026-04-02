#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  check-goodruns-campaigns.sh

Checks the campaign SAM definitions hardcoded in this script and classifies each
definition as:

  YES     explicit good-runs filtering is present in the definition dimensions
  VERIFY  only "beam_good" appears; inspect the full dimensions line manually
  NO      no explicit good-runs filtering was found
  ERROR   samweb describe-definition failed

The script looks for:
  - goodruns: 1
  - defname:goodruns_mcc9...

Notes:
  - Definition names ending in *_goodruns or *_good_runs are hints, not proof.
  - "beam_good" is treated as VERIFY rather than YES.
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v samweb >/dev/null 2>&1; then
  echo "Error: samweb is not available in PATH." >&2
  exit 1
fi

datasets=(
  # Run 1A
  run1_beamon_opentrigger_numi_pandora_3v_a_reprocessing_v08_00_00_84_reco2_reco2_beam_good_good_runs

  # Run 1 beam data
  prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2
  prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2

  # Run 1 FHC
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_all_snapshot
  prodgenie__numi_nue_overlay_mcc9_v08_00_00_48_CV_reco2_run1_reco2
  prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot
  run1_numi_fhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_reco2_reco2

  # Run 1 detvar nu overlay
  prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2

  # Run 1 detvar nue overlay
  prodgenie_numi_intrinsic_nue_Reco2_LYRayleigh_v08_00_00_48_rerun_reco2_run1_reco2
  NuMI_run1_intrinsic_nue_reco2_LYDown_reco2_LYDown_reco2
  prodgenie_numi_overlay_mcc9_v08_00_00_48_TPCDetVar_run1_reco2_nue_SCE_reco2
  prodgenie_numi_intrinsic_nue_Reco2_run1_LYAttenuation_v08_00_00_48_reco2_run1_reco2
  prodgenie_numi_overlay_mcc9_v08_00_00_48_TPCDetVar_run1_reco2_nue_Recomb2_reco2
  prodgenie_numi_nue_Reco2_run1_WireModX_v08_00_00_48_fixedinput_reco2_run1_reco2
  prodgenie_numi_nue_intrinsic_Reco2_run1_WireModYZ_v08_00_00_48_fixed_input_reco2_run1_reco2
  prodgenie__numi_nue_Reco2_run1_WireModThetaXZ_v08_00_00_48_reco2_reco2_run1_reco2
  prodgenie__numi_nue_Reco2_run1_WireModThetaYZ_with_sigma_splines_v08_00_00_48_reco2_reco2_run1_reco2

  # Run 1 RHC
  prodgenie_numi_rhc_nu_overlay_v08_00_00_54_run1_reco2_reco2
  prodgenie_numi_rhc_intrinsic_nue_overlay_v08_00_00_54_run1_reco2_reco2_reco2
  prod_extunbiased_numi_rhc_dirt_overlay_run1_reco2_v08_00_00_67_reco2

  # Run 2a FHC
  prodgenie_numi_overlay_detvar_CV_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2
  run2_numi_nu_overlay_pandora_unified_reco2_run2a_fhc_reco2
  prodgenie_numi_nue_overlay_detvar_CV_run2_FHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_fhc_dirt_overlay_pandora_reco2_run2_reco2
  prod_extnumi_swizzle_inclusive_v4_run2a_reco2_run2a_all_reco2
  prod_numi_swizzle_inclusive_v4_run2_reco2_run2a_beam_good_reco2

  # Run 2a detvar nu overlay
  prodgenie_numi_overlay_detvar_LYAttenuation_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_overlay_detvar_LYRayleigh_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_overlay_detvar_LYDown_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_SCE_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_overlay_detvar_Recomb2_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_detvar_WireModX_standard_FHC_nu_overlay_400k_run2_reco2_reco2_reco2
  prodgenie_numi_overlay_detvar_WireModYZ_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_pandora_reco2_run2_FHC_WireModThetaYZ_withSplines_reco2

  # Run 2a detvar nue overlay
  prodgenie_numi_nue_overlay_detvar_LYAttenuation_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nue_overlay_detvar_LYRayleigh_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nue_overlay_detvar_LYDown_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nue_overlay_detvar_SCE_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nue_overlay_pandora_reco2_run2_FHC_Recomb2_reco2
  prodgenie_numi_nue_overlay_detvar_WireModX_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_overlay_detvar_CV_run2_FHC_intrinsic_nue_WireModYZ_reco2_run2_reco2
  prodgenie_numi_nue_overlay_detvar_WireModThetaXZ_run2_FHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nue_overlay_pandora_reco2_run2_FHC_WireModThetaYZ_withSplines_reco2

  # Run 2b RHC
  prodgenie_numi_nu_overlay_detvar_CV_run2_RHC_reco2_v08_00_00_55_run2_reco2
  run2_numi_nu_overlay_pandora_unified_reco2_run2b_rhc_reco2
  prodgenie_numi_overlay_detvar_CV_run2_RHC_nue_reco2_v08_00_00_55_run2_run2_reco2_goodruns
  prodgenie_numi_rhc_dirt_overlay_pandora_reco2_run2_reco2
  prod_extnumi_swizzle_crt_inclusive_v4b_offbeam_run2_reco2_new2_run2_reco2_all
  prod_numi_swizzle_inclusive_v4b_run2_beamon_run2_reco2_run2_beam_good_reco2

  # Run 2b detvar nu overlay
  prodgenie_numi_nu_overlay_detvar_LYAttenuation_run2_RHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYRayleigh_run2_reco2
  prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYDown_run2_reco2
  prodgenie_numi_overlay_run2_RHC_standard_nu_reco2_detvar_SCE_run2_reco2
  prodgenie_numi_nu_overlay_detvar_Recomb2_run2_RHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_WireModX_run2_RHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_WireModYZ_run2_RHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_RHC_reco2_v08_00_00_55_run2_reco2
  prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_wsplines_run2_RHC_reco2_v08_00_00_5_run2_reco2

  # Run 2b detvar nue overlay
  prodgenie_numi_overlay_detvar_LYAttenuation_run2_RHC_nue_reco2_v08_00_00_55_run2_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_LYRayleigh_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_LYDown_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_SCE_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_Recomb2_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_WireModX_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_WireModYZ_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_detvar_WireModThetaXZ_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns
  prodgenie_numi_nue_overlay_pandora_reco2_run2_RHC_WireModThetaYZ_withSplines_reco2_goodruns

  # Run 3b RHC
  prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_all_snapshot
  prodgenie_numi_intrinsic_nue_overlay_Reco2_CV_run3b_v08_00_00_48_reco2_run3b_reco2
  prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2
  prod_numi_mcc9_v08_00_00_45_run3b_run3b_reco2_beam_good_reco2
  run3_numi_rhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_fix_reco2_reco2

  # Run 3b detvar nu overlay
  numi_run3_standard_nu_overlay_cv_reco2_v08_00_00_54_run3_reco2
  numi_run3_standard_nu_overlay_LYattenuation_reco2_run3_reco2
  numi_run3_standard_nu_overlay_LYrayleigh_reco2_run3_reco2
  numi_run3_standard_nu_overlay_LYDown_reco2_v08_00_00_54_run3_reco2
  numi_run3_standard_nu_overlay_SCE_reco2_v08_00_00_54_run3_reco2
  numi_run3_standard_nu_overlay_Recomb2_reco2_v08_00_00_54_run3_reco2
  numi_run3_standard_nu_overlay_cv_v08_00_00_54_run3_reco2_WiremodX_run3b_reco2
  numi_run3_standard_nu_overlay_WireModYZ_reco2_run3b_reco2
  numi_run3_standard_nu_overlay_WireModThetaXZ_reco2_v08_00_00_54_run3_reco2
  numi_run3_standard_nu_overlay_WireModThetaYZwithsplines_reco2_run3_reco2

  # Run 3b detvar nue overlay
  prodgenie_numi_intrinsic_nue_mcc9_v08_00_00_48_LYRayleigh_run3b_run3b_reco2
  prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_LYDown_reco2_run3b_reco2
  prodgenie_numi_run3b_reco2_LY_Attenuation_intrinsic_nue_run_3b_reco2
  prodgenie_numi_nue_overlay_Reco2_SCE_run3b_v08_00_00_48_reco2_run3b_reco2
  prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_Recomb2_reco1_run3b_reco2_run3b_reco2
  prodgenie_numi_standard_intrinsic_nue_mcc9_v08_00_00_48_WireModX_reco2_run3b_run3b_reco2
  prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_WireModYZ_reco2_run3b_reco2
  prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_WireModThetaXZ_reco2_run3b_reco2
  prodgenie_numi_run3b_reco2_WireModThetaYZ_sigmaSplines_fixed2_intrinsic_nue_run_3b_reco2

  # Run 4a
  run4a_NuMI_RHC_nu_overlay_pandora_unified_reco2_reco2_reco2
  numi_run4a_rhc_nue_overlay_pandora_unified_reco2_run4a_rhc_reco2
  run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_reco2
  NuMI_beam_off_reco_2_for_run_4a_run4a_reco2

  # Run 4b
  numi_run4b_rhc_nu_overlay_pandora_unified_reco2_run4b_reco2
  numi_run4b_rhc_nue_overlay_pandora_unified_reco2_run4b_rhc_reco2
  numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_reco2
  Run_4b_NuMI_beam_off_Pandora_Reco2_run4b_reco2_all
  Run_4b_NuMI_beam_on_Pandora_Reco2_run4b_reco2_beam_good

  # Run 4c
  numi_run4c_fhc_nu_overlay_pandora_unified_reco2_run4c_reco2
  numi_run4c_fhc_nue_overlay_pandora_unified_reco2_run4c_fhc_reco2
  numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_reco2
  run4c_numi_beam_off_pandora_reco2_run4c_reco2_all
  Run_4c_NuMI_beam_on_Pandora_Reco2_run4c_reco2_beam_good

  # Run 4d
  run4d_numi_fhc_nu_overlay_pandora_unified_reco2_run4d_fhc_reco2
  run4d_numi_fhc_nue_overlay_pandora_unified_reco2_run4d_fhc_reco2
  numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_reco2
  run4d_numi_beam_off_pandora_reco2_run4d_reco2_all
  run4d_NuMI_beam_on_pandora_unified_reco2_run4d_reco2_beam_good

  # Run 5
  NuMI_Run5_FHC_nu_overlay_pandora_Reco2_run5_reco2_nonzerolifetime_goodruns
  run5_numi_fhc_nue_overlay_pandora_unified_reco2_run5_fhc_reco2_goodruns
  numi_fhc_dirt_overlay_pandora_reco2_run5_run5_fhc_reco2
  run5_numi_fhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_fix_reco2_reco2
  run5_numi_beam_off_pandora_reco2_all_run5_reco2_all
  run5_numi_beam_on_pandora_reco2_all_run5_reco2_beam_good

  # Strangeness
  prod_strange_resample_fhc_run1_fhc_reco2_reco2_goodruns
  prod_strange_resample_fhc_run2_fhc_reco2_reco2_goodruns
  prod_strange_resample_rhc_run2_rhc_reco2_reco2_goodruns
  prod_strange_resample_rhc_run3_rhc_reco2_reco2_goodruns
)

yes_count=0
verify_count=0
no_count=0
error_count=0

extract_dimensions() {
  awk '
    BEGIN {
      IGNORECASE = 1
      capture = 0
    }
    capture {
      if ($0 ~ /^[[:alpha:]_][[:alnum:]_ -]*[[:space:]]*[:=][[:space:]]*/) exit
      if ($0 ~ /^[[:space:]]*$/) exit
      print
      next
    }
    /^(dimensions|dims)[[:space:]]*[:=][[:space:]]*/ {
      sub(/^[^:=]*[:=][[:space:]]*/, "", $0)
      print
      capture = 1
    }
  ' | paste -sd' ' - | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//'
}

for def in "${datasets[@]}"; do
  desc="$(samweb describe-definition "$def" 2>&1)" || {
    printf 'ERROR   %s\n' "$def"
    printf '  %s\n\n' "$desc"
    error_count=$((error_count + 1))
    continue
  }

  dims="$(printf '%s\n' "$desc" | extract_dimensions)"
  search_blob="$desc"
  note=""
  if [[ -n "$dims" ]]; then
    search_blob="$dims"
  else
    note="(raw describe-definition format did not expose a parsable dimensions line)"
  fi

  if printf '%s\n' "$search_blob" | grep -Eqi 'goodruns:[[:space:]]*1|goodruns_mcc9'; then
    printf 'YES     %s\n\n' "$def"
    yes_count=$((yes_count + 1))
  elif printf '%s\n' "$search_blob" | grep -Eqi 'beam_good'; then
    printf 'VERIFY  %s\n' "$def"
    if [[ -n "$dims" ]]; then
      printf '  %s\n' "$dims"
    elif [[ -n "$note" ]]; then
      printf '  %s\n' "$note"
    fi
    printf '\n'
    verify_count=$((verify_count + 1))
  else
    printf 'NO      %s\n' "$def"
    if [[ -n "$dims" ]]; then
      printf '  %s\n' "$dims"
    elif [[ -n "$note" ]]; then
      printf '  %s\n' "$note"
    fi
    printf '\n'
    no_count=$((no_count + 1))
  fi
done

printf 'Summary: YES=%d VERIFY=%d NO=%d ERROR=%d\n' \
  "$yes_count" "$verify_count" "$no_count" "$error_count"
