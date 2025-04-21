#!/bin/bash

usage() {
    echo "Usage: $0 -r RELEASE -t FILE_TYPE -u USER -n NAME -i INPUTDEF -f FCL -j NUMJOBS -o OUTPUT_XML"
    echo "Options:"
    echo "  -r  Release version"
    echo "  -t  File type"
    echo "  -u  Username"
    echo "  -n  Project name"
    echo "  -i  Input dataset definition"
    echo "  -f  FCL file"
    echo "  -j  Number of jobs"
    echo "  -o  Output XML file"
    exit 1
}

while getopts "r:t:u:n:i:f:j:o:h" opt; do
    case $opt in
        r) RELEASE=$OPTARG ;;
        t) FILE_TYPE=$OPTARG ;;
        u) USER=$OPTARG ;;
        n) NAME=$OPTARG ;;
        i) INPUTDEF=$OPTARG ;;
        f) FCL=$OPTARG ;;
        j) NUMJOBS=$OPTARG ;;
        o) OUTPUT_XML=$OPTARG ;;
        h) usage ;;
        ?) usage ;;
    esac
done

[ -z "$RELEASE" ] || [ -z "$FILE_TYPE" ] || [ -z "$USER" ] || [ -z "$NAME" ] || \
[ -z "$INPUTDEF" ] || [ -z "$FCL" ] || [ -z "$NUMJOBS" ] || [ -z "$OUTPUT_XML" ] && usage

TEMPLATE="../templates/analysis.xml"
[ ! -f "$TEMPLATE" ] && { echo "Error: Template $TEMPLATE not found"; exit 1; }

sed -e "s|{{release}}|$RELEASE|g" \
    -e "s|{{file_type}}|$FILE_TYPE|g" \
    -e "s|{{user}}|$USER|g" \
    -e "s|{{name}}|$NAME|g" \
    -e "s|{{inputdef}}|$INPUTDEF|g" \
    -e "s|{{fcl}}|$FCL|g" \
    -e "s|{{numjobs}}|$NUMJOBS|g" \
    "$TEMPLATE" > "$OUTPUT_XML" || { echo "Error: XML generation failed"; exit 1; }

echo "Generated: $OUTPUT_XML"