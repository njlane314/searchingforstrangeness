#!/bin/bash

set -e

# Default installation directory is $MRB_INSTALL. Allow override via -d.
dir="${MRB_INSTALL}"
while getopts "d:" opt; do
  case "$opt" in
    d)
      dir="$OPTARG"
      ;;
  esac
done
shift $((OPTIND-1))

if [ -z "$dir" ]; then
  echo "Installation directory not specified. Set MRB_INSTALL or use -d." >&2
  exit 1
fi

mkdir -p "$dir"

cp run_strangeness_inference.sh run_inference.py binary_classifier_resnet34.pth "$dir"

ls "$dir"

make_tar_uboone.sh -d "$dir" strangeness.tar
cp -f strangeness.tar /pnfs/uboone/resilient/users/nlane/tarballs/strangeness.tar
rm strangeness.tar

