#!/bin/bash
#
# This script submits workflow stages for beam, or strange.
# It can submit 'reweight' stages, or 'selection' stages.
#
# Usage:
#   ./your_script_name.sh [sample_type] [stage_type]
#
#   sample_type: 'beam' | 'strange'
#   stage_type:  'reweight' | 'selection'
#
# Examples:
#   ./your_script_name.sh strange selection
#   ./your_script_name.sh beam reweight

WORKFLOW_XML="numi_fhc_workflow.xml"

submit_stage() {
    STAGE_NAME=$1
    echo "Processing stage: $STAGE_NAME"
    project.py --xml "${WORKFLOW_XML}" --stage "${STAGE_NAME}" --clean
    project.py --xml "${WORKFLOW_XML}" --stage "${STAGE_NAME}" --submit
    echo "Waiting 1 seconds before next submission..."
    sleep 1
}

check_previous_stage() {
    PREVIOUS_STAGE=$1
    echo "--- Checking status of previous stage: $PREVIOUS_STAGE ---"
    if ! project.py --xml "${WORKFLOW_XML}" --stage "${PREVIOUS_STAGE}" --check; then
        echo "WARNING: Check for stage '$PREVIOUS_STAGE' returned a non-zero exit code."
        echo "Please review its status. Continuing submission, but consider pausing."
    fi
    echo "--- Check complete for $PREVIOUS_STAGE. Review output above. ---"
}

REWEIGHT_BEAM_STAGES=(
    #"reweight_numi_fhc_run1_beam"
)

REWEIGHT_STRANGE_STAGES=(
    #"reweight_numi_fhc_run1_strangeness"
)

SELECTION_BEAM_STAGES=(
    #"selection_numi_fhc_run1_beam_TEST"
    "selection_numi_fhc_run1_beam"
    #"selection_detvar_cv"
    #"selection_detvar_ly_suppression75attenuation8m"
    #"selection_detvar_ly_rayleigh"
    #"selection_detvar_lydown"
    #"selection_detvar_sce"
    #"selection_detvar_recomb2"
    #"selection_detvar_wiremodx"
    #"selection_detvar_wiremodyz"
    #"selection_detvar_wiremodthetaxz"
    #"selection_detvar_wiremodthetayz_withsplines"
)

SELECTION_STRANGE_STAGES=(
    "selection_numi_fhc_run1_strangeness"
    #"selection_detvar_cv_strangeness"
    #"selection_detvar_sce_strangeness"
    #"selection_detvar_recomb2_strangeness"
    #"selection_detvar_ly_down_strangeness"
    #"selection_detvar_ly_rayleigh_strangeness"
    #"selection_detvar_wiremodx_sce_strangeness"
    #"selection_detvar_wiremodyz_sce_strangeness"
    #"selection_detvar_wiremod_yz_strangeness"
    #"selection_detvar_wiremod_thetayz_strangeness" 
)

declare -A SELECTION_INPUT_STAGE_MAP
SELECTION_INPUT_STAGE_MAP["selection_numi_fhc_run1_beam"]="reweight_numi_fhc_run1_beam"
SELECTION_INPUT_STAGE_MAP["selection_numi_fhc_run1_strangeness"]="reweight_numi_fhc_run1_strangeness"

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
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', or 'selection'."
                exit 1
                ;;
        esac
        ;;
    "beam")
        echo "Preparing to submit BEAM stages..."
        case "$STAGE_TYPE" in
            "reweight")
                STAGES_TO_SUBMIT=("${REWEIGHT_BEAM_STAGES[@]}")
                ;;
            "selection")
                STAGES_TO_SUBMIT=("${SELECTION_BEAM_STAGES[@]}")
                ;;
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', or 'selection'."
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
            *)
                echo "Unknown stage type: $STAGE_TYPE. Use 'reweight', or 'selection'."
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
    if [[ "$STAGE_TYPE" == "selection" ]]; then
        PREVIOUS_REWEIGHT_STAGE="${SELECTION_INPUT_STAGE_MAP[$stage]}"
        if [[ -n "$PREVIOUS_REWEIGHT_STAGE" ]]; then
            echo "Identified '$stage' as dependent on '$PREVIOUS_REWEIGHT_STAGE'."
            check_previous_stage "$PREVIOUS_REWEIGHT_STAGE"
        else
            echo "Stage '$stage' is a selection stage but does not have a specific reweight input stage defined in the map (e.g., it might use inputdef)."
        fi
    fi
    submit_stage "${stage}"
done

echo "All selected stages submitted."