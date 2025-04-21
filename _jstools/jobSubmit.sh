#!/bin/bash

usage() {
  echo "Usage: $0 XML_FILE"
  echo "  XML_FILE: Path to the XML file to submit"
  exit 1
}

[ $# -ne 1 ] && usage
XML_FILE="$1"

[ ! -f "$XML_FILE" ] && { echo "Error: $XML_FILE not found"; exit 1; }

echo "Submitting $XML_FILE..."
project.py --xml "$XML_FILE" --stage analyse --submit > submission.log 2>&1
[ $? -eq 0 ] && echo "Submitted successfully. See submission.log" || { echo "Submission failed"; exit 1; }