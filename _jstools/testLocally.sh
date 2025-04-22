#!/bin/bash

usage() {
    echo "Usage: runFhiclTest [OPTIONS] FICLFILE"
    echo
    echo "Retrieve a ROOT file from dCache and test if a FHiCL configuration runs successfully."
    echo
    echo "Arguments:"
    echo "  FICLFILE  The FHiCL configuration file to test (e.g., selectionconfig.fcl)"
    echo
    echo "Options:"
    echo "  -s SAMDEF  Specify the SAM definition to use (default: prod_strange_resample_fhc_run2_fhc_reco2_reco2)"
    echo "  -o OUTDIR  Specify the output directory (default: /exp/uboone/data/users/$USER/analysis)"
    echo "  -v         Enable verbose output"
    echo "  -h         Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - Valid x509 proxy (run 'kinit' and 'voms-proxy-init -voms uboone:/uboone/Role=Analysis')"
    echo "  - samweb and lar commands available in the environment"
    echo
    echo "Examples:"
    echo "  runFhiclTest selectionconfig.fcl"
    echo "  runFhiclTest -s my_dataset -v -o /path/to/output selectionconfig.fcl"
}

SAMDEF="prod_strange_resample_fhc_run2_fhc_reco2_reco2"
OUTDIR="/exp/uboone/data/users/$USER/analysis"
VERBOSE=""

while getopts s:o:vh OPT; do
    case $OPT in
        s)
            SAMDEF="$OPTARG"
            ;;
        o)
            OUTDIR="$OPTARG"
            ;;
        v)
            VERBOSE=yes
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            echo "unknown option, exiting"
            usage
            exit 1
            ;;
    esac
done

shift $(expr $OPTIND - 1)
FICLFILE="$1"

if [ -z "$FICLFILE" ]; then
    echo "ERROR - one FHiCL file argument is required"
    usage
    exit 1
fi

if [ ! -f "$FICLFILE" ]; then
    echo "ERROR - FHiCL file '$FICLFILE' does not exist"
    exit 1
fi

FINV=/tmp/x509up_u$UID
[ -n "$X509_USER_PROXY" ] && FINV="$X509_USER_PROXY"

if [ ! -e "$FINV" ]; then
    echo "ERROR - could not find x509 proxy. Run 'kinit' and 'voms-proxy-init -voms uboone:/uboone/Role=Analysis'"
    exit 1
fi

if ! command -v samweb > /dev/null 2>&1; then
    echo "ERROR - samweb not found. Please set up the MicroBooNE environment."
    exit 1
fi

if ! command -v lar > /dev/null 2>&1; then
    echo "ERROR - lar command not found. Please set up the art framework environment."
    exit 1
fi

if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
    if [ $? -ne 0 ]; then
        echo "ERROR - could not create output directory '$OUTDIR'"
        exit 1
    fi
fi

echo "Fetching a single file from SAM definition: $SAMDEF"
FILE=$(samweb list-files "defname:$SAMDEF" | head -n 1)

if [ -z "$FILE" ]; then
    echo "ERROR - no files found in SAM definition '$SAMDEF'"
    exit 1
fi

[ "$VERBOSE" ] && echo "Selected file: $FILE"

LOCATION=$(samweb locate-file "$FILE" | grep -o '/pnfs/.*' | head -n 1)

if [ -z "$LOCATION" ]; then
    echo "ERROR - could not locate dCache path for file '$FILE'"
    exit 1
fi

FILEPATH="${LOCATION}/${FILE}"

if [ ! -f "$FILEPATH" ]; then
    echo "ERROR - file not found at '$FILEPATH'"
    exit 1
fi

[ "$VERBOSE" ] && echo "File path: $FILEPATH"

OUTPUTFILE="${OUTDIR}/test_output.root"
LOGFILE="${OUTDIR}/lar_log.txt"

echo "Running: lar -c $FICLFILE -s $FILEPATH -T $OUTPUTFILE"
if [ "$VERBOSE" ]; then
    lar -c "$FICLFILE" -s "$FILEPATH" -T "$OUTPUTFILE" 2>&1 | tee "$LOGFILE"
    STATUS=${PIPESTATUS[0]}
else
    lar -c "$FICLFILE" -s "$FILEPATH" -T "$OUTPUTFILE" > "$LOGFILE" 2>&1
    STATUS=$?
fi

if [ $STATUS -ne 0 ]; then
    echo "ERROR - FHiCL processing failed. See $LOGFILE for details."
    exit 1
elif [ ! -f "$OUTPUTFILE" ]; then
    echo "ERROR - output file '$OUTPUTFILE' was not created."
    exit 1
else
    echo "Success - FHiCL configuration ran successfully. Output: $OUTPUTFILE"
fi

exit 0