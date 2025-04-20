#!/bin/bash

usage() {
    echo "Usage: samPrestageWorkflow [OPTIONS]"
    echo
    echo "Prestage multiple SAM definitions from tape to dCache in MicroBooNE."
    echo "Reads definitions from a file (default: definitions.txt)."
    echo
    echo "Options:"
    echo "  -f FILE  File with SAM definitions (one per line)"
    echo "  -v       Enable verbose output"
    echo "  -h       Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - MicroBooNE environment"
    echo "  - sam-web-client"
    echo "  - Valid x509 proxy"
    echo
    echo "Examples:"
    echo "  samPrestageWorkflow"
    echo "  samPrestageWorkflow -f my_defs.txt -v"
}

prestage_definition() {
    local DEF="$1"
    local VERBOSE="$2"
    echo "Starting prestage for '$DEF' at $(date)"

    if ! samweb describe-definition "$DEF" > /dev/null 2>&1; then
        echo "Warning: Definition '$DEF' does not exist. Skipping."
        return 1
    fi

    local SAM_PROJECT="prestage_${DEF}_$(date +%s)"
    local SAM_PROJECT_URL=$(samweb start-project --defname="$DEF" "$SAM_PROJECT")
    if [ -z "$SAM_PROJECT_URL" ]; then
        echo "Error: Failed to start SAM project for '$DEF'."
        return 1
    fi
    [ "$VERBOSE" ] && echo "SAM project started: $SAM_PROJECT"

    local SAM_CONSUMER_ID=$(samweb start-process --appname=null --appversion=0 "$SAM_PROJECT_URL")
    if [ -z "$SAM_CONSUMER_ID" ]; then
        echo "Error: Failed to start consumer for '$DEF'."
        return 1
    fi
    [ "$VERBOSE" ] && echo "Consumer process started: $SAM_CONSUMER_ID"

    local N=0
    local SAM_FILE=$(samweb get-next-file "$SAM_PROJECT_URL" "$SAM_CONSUMER_ID")
    while [ -n "$SAM_FILE" ]; do
        N=$((N + 1))
        samweb release-file --status=ok "$SAM_PROJECT_URL" "$SAM_CONSUMER_ID" "$SAM_FILE"
        (( N % 100 == 0 )) && echo "Processed $N files for '$DEF'"
        SAM_FILE=$(samweb get-next-file "$SAM_PROJECT_URL" "$SAM_CONSUMER_ID")
    done

    echo "Completed prestaging $N files for '$DEF' at $(date)"
}

DEFS_FILE="definitions.txt"
VERBOSE=""

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

while [ "$1" ]; do
    if [ "$1" == "-f" ]; then
        shift
        DEFS_FILE="$1"
    elif [ "$1" == "-v" ]; then
        VERBOSE="yes"
    else
        echo "Unknown option: $1"
        usage
        exit 1
    fi
    shift
done

if [ ! -f "$DEFS_FILE" ]; then
    echo "Error: Definitions file '$DEFS_FILE' not found."
    exit 1
fi

if ! command -v samweb &> /dev/null; then
    echo "Please set up sam-web-client (e.g., via uboone_setup)"
    exit 1
fi

FINV="/tmp/x509up_u$UID"
[ -n "$X509_USER_PROXY" ] && FINV="$X509_USER_PROXY"
if [ ! -e "$FINV" ]; then
    echo "Please set up your x509 proxy (e.g., with 'kinit' and 'voms-proxy-init')"
    exit 1
fi

echo "Reading definitions from '$DEFS_FILE'"
cat "$DEFS_FILE" | while read DEF; do
    [ -z "$DEF" ] && continue

    echo "----------------------------------------"
    echo "Processing definition: '$DEF'"
    prestage_definition "$DEF" "$VERBOSE" || echo "Failed to prestage '$DEF'. Moving to next."
    echo "----------------------------------------"
done

echo "Prestage workflow completed at $(date)"
exit 0