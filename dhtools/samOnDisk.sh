#!/bin/bash

usage() {
    echo "Usage: samOnDisk [OPTIONS] DATASET"
    echo
    echo "Check how many files in a MicroBooNE dataset are on disk in dCache."
    echo "Stops after 1000 checks or total files, whichever is smaller."
    echo
    echo "Arguments:"
    echo "  DATASET  Dataset name or definition"
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
    echo "  samOnDisk -v my_dataset"
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

TMP=$(mktemp)
echo "Listing files..."
samweb list-files "dataset_def_name=$DS" > $TMP
NT=$(wc -l < $TMP)
echo "Found $NT files"

NN=1
ND=0
while [[ $NN -le 1000 && $NN -le $NT ]]; do
    if [ $NT -lt 1000 ]; then
        SF=$(awk "NR==$NN" $TMP)
    else
        SF=$(shuf -n 1 $TMP)
    fi
    SLFS=$(samweb locate-file "$SF" | grep enstore)
    PNFSDIR=$(echo "$SLFS" | awk -F: '{print $2}' | sed 's/(.*)//')
    ANS=$(cat "$PNFSDIR/'.(get)('$SF')(locality)'")
    if [[ "$ANS" =~ "ONLINE" ]]; then
        ND=$((ND + 1))
    fi
    RR=$(echo "$ND $NN" | awk '{print $1*100/$2}')
    printf "%4d/%4d are on disk, %4.1f %%\n" $ND $NN $RR
    NN=$((NN + 1))
done

rm -f $TMP
exit 0