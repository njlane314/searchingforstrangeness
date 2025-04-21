#!/bin/bash

usage() {
    echo "Usage: samListLocations [OPTIONS] SAMWEB-ARGS..."
    echo
    echo "List file locations for a SAM dataset in MicroBooNE, reformatted as dCache filespecs."
    echo
    echo "Arguments:"
    echo "  SAMWEB-ARGS  Arguments for 'samweb list-file-locations' (e.g., --defname=NAME)"
    echo
    echo "Options:"
    echo "  -d        Print only files on disk (requires x509 proxy)"
    echo "  -f        Print files on disk first, then others (requires x509 proxy)"
    echo "  -v        Enable verbose output"
    echo "  -h        Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - sam-web-client"
    echo "  - For -d or -f: valid x509 proxy"
    echo
    echo "Examples:"
    echo "  samListLocations --defname=my_dataset"
    echo "  samListLocations -d --dim \"dh.dataset=my_dataset\""
}

get_locality() {
    local filespec="$1"
    local cfilespec=$(echo "$filespec" | sed 's|/pnfs|/pnfs/fnal.gov/usr|')
    local finv=/tmp/x509up_u$UID
    [ -n "$X509_USER_PROXY" ] && finv="$X509_USER_PROXY"
    local ans=$(curl -s -L --capath /etc/grid-security/certificates --cert "$finv" --cacert "$finv" --key "$finv" -X GET "https://fndca1.fnal.gov:3880/api/v1/namespace${cfilespec}?locality=true")
    local loc=$(echo "$ans" | grep fileLocality | tr -d '":,' | awk '{print $2}')
    echo "$loc"
}

DODISK=""
DOFIRST=""
DOLOC=""

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

TMPF=""
if [ "$1" == "-d" ]; then
    shift
    DODISK=yes
    DOLOC=yes
elif [ "$1" == "-f" ]; then
    shift
    DOFIRST=yes
    DOLOC=yes
    TMPF=$(mktemp)
fi

if ! command -v samweb &> /dev/null; then
    echo "Please set up sam-web-client (e.g., via uboone_setup)"
    exit 1
fi

if [ "$DOLOC" ]; then
    FINV=/tmp/x509up_u$UID
    [ -n "$X509_USER_PROXY" ] && FINV="$X509_USER_PROXY"
    if [ ! -e "$FINV" ]; then
        echo "Please set up your x509 proxy (e.g., with 'kinit' and 'voms-proxy-init')"
        exit 1
    fi
fi

TMPL=$(mktemp)

samweb list-file-locations "${@}" | \
awk '{if(NF==4) print $2" "$3" "$4" "$1; else print $1" "$2" "$3;}' | \
sed -e 's/enstore:/enstore /' -e 's/dcache:/dcache /' > $TMPL

cat $TMPL | awk '{print $3}' | sort | uniq | while read FILENAME; do
    FSS=$(cat $TMPL | awk '{if($3=="'$FILENAME'") print $1" "$2" "$5}' | sort)

    TAPESPEC=$(echo "$FSS" | awk '{ if($1=="enstore") print $2"/'$FILENAME'" }')
    DISKSPEC=$(echo "$FSS" | awk '{ if($1=="dcache") print $2"/'$FILENAME'" }')

    TAPEURI=$(echo "$FSS" | awk '{ if($1=="enstore") print $3 }')
    [ -z "$TAPEURI" ] && TAPEURI="$TAPESPEC"
    DISKURI=$(echo "$FSS" | awk '{ if($1=="dcache") print $3 }')
    [ -z "$DISKURI" ] && DISKURI="$DISKSPEC"

    FILEURI="none"
    if [ "$DISKSPEC" ]; then
        FILESPEC="$DISKSPEC"
        FILEURI="$DISKURI"
        LOC="ONLINE"
    elif [ "$TAPESPEC" ]; then
        FILESPEC="$TAPESPEC"
        FILEURI="$TAPEURI"
        LOC=""
        if [ "$DOLOC" ]; then
            LOC=$(get_locality "$FILESPEC")
        fi
    fi

    if [ "$DODISK" ]; then
        [[ "$LOC" == "ONLINE" ]] && echo "$FILEURI"
    elif [ "$DOFIRST" ]; then
        if [[ "$LOC" == "ONLINE" ]]; then
            echo "$FILEURI"
        else
            echo "$FILEURI" >> "$TMPF"
        fi
    else
        echo "$FILEURI"
    fi
done

if [ "$DOFIRST" ]; then
    cat "$TMPF"
    rm -f "$TMPF"
fi
rm -f "$TMPL"

exit 0