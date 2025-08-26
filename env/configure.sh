#!/bin/bash
if [[ -z "${STRANGENESS_DIR}" ]]; then
  echo "STRANGENESS_DIR is not set. Source env/setenv.sh first."
  return 1 2>/dev/null || exit 1
fi

source /usr/local/root/bin/thisroot.sh

