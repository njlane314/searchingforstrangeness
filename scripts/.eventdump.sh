#!/bin/sh

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $(basename "$0") <sam_definition>"
    exit 1
fi

SAM_DEF="$1"
FHICL_FILE="eventdump.fcl"

echo "Fetching one file from SAM definition: ${SAM_DEF}"
file=$(samweb list-files defname:"${SAM_DEF}" | head -n 1)

if [ -z "${file}" ]; then
    echo "Error: No files found for SAM definition '${SAM_DEF}'."
    exit 1
fi

echo "Locating file: ${file}"
filedir=$(samweb locate-file "${file}" | grep -o '/pnfs/.*' | sed 's/(\([0-9]*@[a-z0-9]*\))//g' | head -n 1)

if [ -z "${filedir}" ]; then
    echo "Error: Could not locate the file directory."
    exit 1
fi

filepath="${filedir}/${file}"

if [ ! -f "${filepath}" ]; then
    echo "Error: File not found at '${filepath}'."
    exit 1
fi

echo "Dumping events from: ${filepath}"
lar -c "${FHICL_FILE}" -s "${filepath}" -n 1

echo "Event dump complete."