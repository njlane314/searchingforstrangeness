#!/bin/bash

usage() {
    echo "Usage: samPrestage -d DEF | -s DS [OPTIONS]"
    echo
    echo "Prestage files from tape to dCache for a MicroBooNE dataset or definition."
    echo
    echo "Options:"
    echo "  -d DEF   Specify dataset definition"
    echo "  -s DS    Specify dataset name"
    echo "  -v       Enable verbose output"
    echo "  -h       Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - MicroBooNE environment"
    echo "  - sam-web-client"
    echo "  - Valid x509 proxy"
    echo
    echo "Examples:"
    echo "  samPrestage -d my_definition"
    echo "  samPrestage -s my_dataset -v"
}

export V=""
export DS=""
export DD=""
while getopts hvd:s: OPT; do
    case $OPT in
        d)
            export DD=$OPTARG
            ;;
        s)
            export DS=$OPTARG
            ;;
        h)
            usage
            exit 0
            ;;
        v)
            export V="true"
            ;;
        *)
            echo unknown option, exiting
            exit 1
            ;;
    esac
done

if [ -z "$DD" ] && [ -z "$DS" ]; then
    usage
    exit 1
fi

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

if [ -n "$DD" ]; then
    export SAM_DD=$DD
else
    export SAM_DD=${USER}_zzz_samPrestage_$(date +%s)
    samweb create-definition $SAM_DD "dataset=$DS"
    samweb take-snapshot $SAM_DD
fi

echo "Summary of files to be prestaged:"
samweb list-files --summary "dataset_def_name=$SAM_DD"

export SAM_PROJECT="prestage_${SAM_DD}_$(date +%s)"
export SAM_PROJECT_URL=$(samweb start-project --defname=$SAM_DD $SAM_PROJECT)

if [ "$V" ]; then
    echo SAM_PROJECT_URL=$SAM_PROJECT_URL
fi
if [ -z "$SAM_PROJECT_URL" ]; then
    echo "SAM_PROJECT_URL not set - exiting!"
    exit 1
fi

export SAM_CONSUMER_ID=$(samweb start-process --appname=null --appversion=0 $SAM_PROJECT_URL)
if [ "$V" ]; then
    echo SAM_CONSUMER_ID=$SAM_CONSUMER_ID
fi
if [ -z "$SAM_CONSUMER_ID" ]; then
    echo "SAM_CONSUMER_ID not set - exiting!"
    exit 1
fi

STIME=$(date +%s)
N=0
export SAM_FILE=$(samweb get-next-file $SAM_PROJECT_URL $SAM_CONSUMER_ID)
while [ -n "$SAM_FILE" ]; do
    N=$((N + 1))
    if [ "$V" ]; then
        echo SAM_FILE=$SAM_FILE
    fi

    samweb release-file --status=ok $SAM_PROJECT_URL $SAM_CONSUMER_ID $SAM_FILE

    if (( N % 10 == 0 && N < 100 )); then
        echo Processing $N
    elif (( N % 100 == 0 && N < 1000 )); then
        echo Processing $N
    elif (( N % 1000 == 0 )); then
        echo Processing $N
    fi

    export SAM_FILE=$(samweb get-next-file $SAM_PROJECT_URL $SAM_CONSUMER_ID)
done

DTIME=$(date +%s)
ELTIME=$((DTIME - STIME))

echo "Prestaged $N files in $ELTIME seconds"

exit 0