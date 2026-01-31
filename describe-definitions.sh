#!/usr/bin/env bash
set -euo pipefail

setup sam_web_client
export SAM_EXPERIMENT=uboone

datasets=(
  # --- Beam ---
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3

  # --- Detector variations (Run 1 FHC) ---
  prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2
  prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2
  prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2

  # --- Strangeness + detvars ---
  prod_strange_resample_fhc_run1_fhc_reco2_reco2
  detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2
  Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2
  Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2
  detvar_prod_strange_resample_fhc_run1_respin_wiremodX_sce_reco2_reco2
  detvar_prod_strange_resample_fhc_run1_respin_wiremodYZ_sce_reco2_reco2
  Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2
  Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2
  detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2
  detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2

  # --- Dirt ---
  prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0
  prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1

  # --- EXT & Data (Run 1 FHC) ---
  prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2
  prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2
)

for def in "${datasets[@]}"; do
  echo "================================================================"
  echo "DATASET: $def"
  echo

  echo ">>> samweb describe-definition $def"
  samweb describe-definition "$def" || { echo "  (describe-definition failed)"; echo; continue; }

  echo
  echo ">>> samweb list-definition-files --summary $def"
  samweb list-definition-files --summary "$def" || echo "  (list-definition-files failed)"

  echo
  echo ">>> samweb list-files --summary \"isparentof:( defname: $def ) and availability: anylocation\""
  samweb list-files --summary "isparentof:( defname: $def ) and availability: anylocation" \
    || echo "  (list-files parents summary failed)"

  echo

done
