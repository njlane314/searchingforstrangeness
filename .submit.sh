#!/bin/bash
#
# This script submits workflow stages for beam, strange, or both sample types.
# It can submit 'reweight' stages, 'selection' stages, or 'all' stages.
#
# Usage:
#   ./your_script_name.sh [sample_type] [stage_type]
#
#   sample_type: 'beam' | 'strange' | (empty for both)
#   stage_type:  'reweight' | 'selection' | 'all' | (empty for all)
#
# Examples:
#   ./your_script_name.sh strange selection  - Submits only strange selection stages.
#   ./your_script_name.sh beam reweight      - Submits only beam reweight stages.
#   ./your_script_name.sh "" all             - Submits all strange and beam reweight/selection stages.
#   ./your_script_name.sh                    - Same as above, "" all is default.
#

WORKFLOW_XML="numi_fhc_workflow.xml"

submit_stage() {
    STAGE_NAME=$1
    echo "Processing stage: $STAGE_NAME"
    project.py --xml "${WORKFLOW_XML}" --stage "${STAGE_NAME}" --clean
    project.py --xml "${WORKFLOW_XML}" --stage "${STAGE_NAME}" --submit
    echo "Waiting 5 seconds before next submission..."
    sleep 5
}

REWEIGHT_BEAM_STAGES=(
    "reweight_numi_fhc_run1_beam"
    "reweight_detvar_cv"
    "reweight_detvar_ly_suppression75attenuation8m"
    "reweight_detvar_ly_rayleigh"
    "reweight_detvar_lydown"
    "reweight_detvar_sce"
    "reweight_detvar_recomb2"
    "reweight_detvar_wiremodx"
    "reweight_detvar_wiremodyz"
    "reweight_detvar_wiremodthetaxz"
    "reweight_detvar_wiremodthetayz_withsplines"
)

REWEIGHT_STRANGE_STAGES=(
    "reweight_numi_fhc_run1_strangeness"
    "reweight_detvar_cv_strangeness"
    "reweight_detvar_ly_rayleigh_strangeness"
    "reweight_detvar_ly_down_strangeness"
    "reweight_detvar_wiremodx_sce_strangeness"
    "reweight_detvar_wiremodyz_sce_strangeness"
    "reweight_detvar_wiremod_yz_strangeness"
    "reweight_detvar_wiremod_thetayz_strangeness"
    "reweight_detvar_sce_strangeness"
    "reweight_detvar_recomb2_strangeness"
)

SELECTION_BEAM_STAGES=(
    "selection_numi_fhc_run1_beam"
    "selection_detvar_cv"
    "selection_detvar_ly_suppression75attenuation8m"
    "selection_detvar_ly_rayleigh"
    "selection_detvar_lydown"
    "selection_detvar_sce"
    "selection_detvar_recomb2"
    "selection_detvar_wiremodx"
    "selection_detvar_wiremodyz"
    "selection_detvar_wiremodthetaxz"
    "selection_detvar_wiremodthetayz_withsplines"
)

SELECTION_STRANGE_STAGES=(
    "selection_numi_fhc_run1_strangeness"
    "selection_detvar_cv_strangeness"
    "selection_detvar_ly_rayleigh_strangeness"
    "selection_detvar_ly_down_strangeness"
    "selection_detvar_wiremodx_sce_strangeness"
    "selection_detvar_wiremodyz_sce_strangeness"
    "selection_detvar_wiremod_yz_strangeness"
    "selection_detvar_wiremod_thetayz_strangeness"
    "selection_detvar_sce_strangeness"
    "selection_detvar_recomb2_strangeness"
)

STAGES_TO_SUBMIT=()
SAMPLE_TYPE="$1"
STAGE_TYPE="$2"

case "$SAMPLE_TYPE" in
    "strange")
        echo "Preparing to submit STRANGENESS stages..."
        case "$STAGE_TYPE" in
            "reweight")
                STAGES_TO_SUBMIT=("${REWEIGHT_STRANGE_STAGES[@]}")
                ;;
            "selection")
                STAGES_TO_SUBMIT=("${SELECTION_STRANGE_STAGES[@]}")
                ;;
            "all" | "")
                STAGES_TO_SUBMIT=("${REWEIGHT_STRANGE_STAGES[@]}" "${SELECTION_STRANGE_STAGES[@]}")
                ;;
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', 'selection', or 'all'."
                exit 1
                ;;
        esac
        ;;
    "beam")
        echo "Preparing to submit BEAM and general DETVAR stages..."
        case "$STAGE_TYPE" in
            "reweight")
                STAGES_TO_SUBMIT=("${REWEIGHT_BEAM_STAGES[@]}")
                ;;
            "selection")
                STAGES_TO_SUBMIT=("${SELECTION_BEAM_STAGES[@]}")
                ;;
            "all" | "")
                STAGES_TO_SUBMIT=("${REWEIGHT_BEAM_STAGES[@]}" "${SELECTION_BEAM_STAGES[@]}")
                ;;
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', 'selection', or 'all'."
                exit 1
                ;;
        esac
        ;;
    "")
        echo "No sample type specified. Preparing to submit both STRANGENESS and BEAM stages."
        case "$STAGE_TYPE" in
            "reweight")
                STAGES_TO_SUBMIT=("${REWEIGHT_STRANGE_STAGES[@]}" "${REWEIGHT_BEAM_STAGES[@]}")
                ;;
            "selection")
                STAGES_TO_SUBMIT=("${SELECTION_STRANGE_STAGES[@]}" "${SELECTION_BEAM_STAGES[@]}")
                ;;
            "all" | "")
                echo "Submitting STRANGENESS reweight stages first..."
                STAGES_TO_SUBMIT=("${REWEIGHT_STRANGE_STAGES[@]}")
                echo ""
                echo "Submitting BEAM and general DETVAR reweight stages now..."
                STAGES_TO_SUBMIT+=("${REWEIGHT_BEAM_STAGES[@]}")
                echo ""
                echo "Submitting STRANGENESS selection stages now..."
                STAGES_TO_SUBMIT+=("${SELECTION_STRANGE_STAGES[@]}")
                echo ""
                echo "Submitting BEAM and general DETVAR selection stages now..."
                STAGES_TO_SUBMIT+=("${SELECTION_BEAM_STAGES[@]}")
                ;;
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', 'selection', or 'all'."
                exit 1
                ;;
        esac
        ;;
    *)
        echo "Unknown sample type: $SAMPLE_TYPE. Use 'beam', 'strange', or leave empty for both."
        exit 1
        ;;
esac

if [ ${#STAGES_TO_SUBMIT[@]} -eq 0 ]; then
    echo "No stages selected for submission. Exiting."
    exit 0
fi

for stage in "${STAGES_TO_SUBMIT[@]}"; do
    submit_stage "${stage}"
done

echo "All selected stages submitted."