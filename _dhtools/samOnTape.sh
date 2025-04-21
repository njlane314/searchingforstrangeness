#!/bin/bash

usage() {
    echo "Usage: samOnTape [OPTIONS] DATASET"
    echo
    echo "Summarize how many files in a MicroBooNE dataset have tape locations."
    echo
    echo "Arguments:"
    echo "  DATASET  Dataset name"
    echo
    echo "Options:"
    echo "  -v       Enable verbose output"
    echo "  -h       Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - MicroBooNE environment (SAM_EXPERIMENT=uboone)"
    echo "  - sam-web-client"
    echo
    echo "Example:"
    echo "  samOnTape -v my_dataset"
}

while getopts h OPT; do
    case $OPT in
        h)
            usage
            exit 0
            ;;
        *)
            echo unknown option, exiting
            usage
            exit 1
            ;;
    esac
done

RC=0
if [ -z "$SAM_EXPERIMENT" ] || [ "$SAM_EXPERIMENT" != "uboone" ]; then
    echo "Please set up the MicroBooNE environment (e.g., source uboone_setup)"
    RC=1
fi
if ! command -v samweb >& /dev/null; then
    echo "Please set up sam-web-client"
    RC=1
fi
[ $RC -ne 0 ] && exit 1

DS="$1"
if [ -z "$DS" ]; then
    echo "ERROR - no dataset provided"
    usage
    exit 1
fi

NT=$(samweb list-files --summary "dataset_def_name=$DS" | grep File | awk '{print $3}')
if [ -z "$NT" ]; then
    echo "ERROR - Failed to retrieve file count for dataset $DS"
    exit 1
fi

NP=$(samweb list-files --summary "dataset_def_name=$DS and tape_label %" | grep File | awk '{print $3}')
if [ -z "$NP" ]; then
    echo "ERROR - Failed to retrieve tape file count for dataset $DS"
    exit 1
fi

NN=$((NT - NP))

printf "%4d Total   %4d on tape  %4d not on tape\n" $NT $NP $NN

exit 0