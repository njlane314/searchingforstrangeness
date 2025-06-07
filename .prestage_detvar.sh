#!/bin/bash

LOG_DIR=".tmp"

PARENT_DEFINITIONS=(
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
)

mkdir -p "${LOG_DIR}"

echo "Starting prestage for the specified SAMWeb definitions."
echo "---"

for PARENT_DEF in "${PARENT_DEFINITIONS[@]}"; do
    samweb describe-definition "${PARENT_DEF}" &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error: Definition '${PARENT_DEF}' does not exist. Skipping."
        echo "---"
        continue
    fi
    
    LOG_FILE="${LOG_DIR}/prestage_${PARENT_DEF}.log"
    
    echo "Starting prestage for definition: ${PARENT_DEF}"
    echo "Logging to: ${LOG_FILE}"
    
    nohup samweb prestage-dataset --defname="${PARENT_DEF}" >& "${LOG_FILE}" &
    
    echo "---"
done

echo "All prestage jobs have been launched in the background."
echo "Monitor progress by checking the log files in the '${LOG_DIR}' directory."