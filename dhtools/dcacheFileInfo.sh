#!/bin/bash

usage() {
    echo "Usage: dcacheFileInfo [OPTIONS] FILESPEC"
    echo
    echo "Retrieve information about a dCache file in the MicroBooNE experiment."
    echo
    echo "Arguments:"
    echo "  FILESPEC  Full dCache filespec (e.g., /pnfs/uboone/path/to/file.root) or filename with -t or -p"
    echo
    echo "Options:"
    echo "  -c        Print dCache checksum"
    echo "  -l        Print locality (NEARLINE, ONLINE, ONLINE_AND_NEARLINE)"
    echo "  -d        Print creation date (epoch seconds)"
    echo "  -t        Interpret FILESPEC as filename and lookup tape location"
    echo "  -p        Interpret FILESPEC as filename and lookup disk location"
    echo "  -v        Enable verbose output"
    echo "  -h        Display this help message"
    echo
    echo "Prerequisites:"
    echo "  - Valid x509 proxy (run 'kinit' and 'voms-proxy-init -voms uboone:/uboone/Role=Analysis')"
    echo
    echo "Examples:"
    echo "  dcacheFileInfo -c -l /pnfs/uboone/path/to/file.root"
    echo "  dcacheFileInfo -t -v filename.root"
}

DOCRC=""
DOLOC=""
DODATE=""
VERBOSE=""
LTAPE=""
LDISK=""

while getopts cldtpvh OPT; do
    case $OPT in
        c)
            DOCRC=yes
            ;;
        l)
            DOLOC=yes
            ;;
        d)
            DODATE=yes
            ;;
        t)
            LTAPE=yes
            ;;
        p)
            LDISK=yes
            ;;
        v)
            VERBOSE=yes
            ;;
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

shift $(expr $OPTIND - 1)
FFS="$@"

if [ -z "$FFS" ]; then
    echo "ERROR - one pnfs filespec or filename argument is required"
    exit 1
fi

FINV=/tmp/x509up_u$UID
[ -n "$X509_USER_PROXY" ] && FINV="$X509_USER_PROXY"

if [ ! -e $FINV ]; then
    echo "ERROR - could not find x509 proxy. Run 'kinit' and 'voms-proxy-init -voms uboone:/uboone/Role=Analysis'"
    exit 1
fi

if [[ -n "$LTAPE" || -n "$LDISK" ]]; then
    if ! command -v samweb > /dev/null 2>&1; then
        echo "ERROR - samweb not found. Please set up the MicroBooNE environment."
        exit 1
    fi
    location=$(samweb locate-file "$FFS")
    if [ -z "$location" ]; then
        echo "ERROR - could not find location for filename $FFS"
        exit 1
    fi
    FFS=$(echo "$location" | sed 's/(.*)//')
    if [ -z "$FFS" ]; then
        echo "ERROR - failed to parse location for $FFS"
        exit 1
    fi
fi

CFFS=$(echo "$FFS" | sed 's|/pnfs|/pnfs/fnal.gov/usr|')

VS=" -s "
[ "$VERBOSE" ] && VS=""

ANS=$(curl $VS -L --capath /etc/grid-security/certificates --cert "$FINV" --cacert "$FINV" --key "$FINV" -X GET "https://fndca1.fnal.gov:3880/api/v1/namespace${CFFS}?checksum=true&locality=true")

RC=$?
if [ $RC -ne 0 ]; then
    echo "ERROR - curl command failed"
    exit 1
fi

[ "$VERBOSE" ] && echo "$ANS"

CRC=""
if [ "$DOCRC" ]; then
    CRC=$(echo -n "$ANS" | grep value | tr -d '":,' | awk '{print $2}')
fi

LOC=""
if [ "$DOLOC" ]; then
    LOC=$(echo -n "$ANS" | grep fileLocality | tr -d '":,' | awk '{print $2}')
fi

DATE=""
if [ "$DODATE" ]; then
    DATEMS=$(echo -n "$ANS" | grep creationTime | tr -d '":,' | awk '{print $2}')
    DATE=$(($DATEMS / 1000))
fi

[ -n "$CRC$LOC$DATE" ] && echo "$CRC $LOC $DATE"

exit 0