#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env/setenv.sh"
source "${SCRIPT_DIR}/env/configure.sh"

# Configure HDF5 include and library paths
if command -v h5c++ >/dev/null 2>&1; then
  HDF5_INC=$(h5c++ -showconfig | awk -F: '/Include directories/ {print $2}' | awk '{print $1}')
  HDF5_LIB=$(h5c++ -showconfig | awk -F: '/Library directories/ {print $2}' | awk '{print $1}')
  export CPLUS_INCLUDE_PATH="${HDF5_INC}:${CPLUS_INCLUDE_PATH}"
  export LIBRARY_PATH="${HDF5_LIB}:${LIBRARY_PATH}"
  export LD_LIBRARY_PATH="${HDF5_LIB}:${LD_LIBRARY_PATH}"
elif [ -d /usr/include/hdf5/serial ] && [ -d /usr/lib/x86_64-linux-gnu/hdf5/serial ]; then
  export CPLUS_INCLUDE_PATH="/usr/include/hdf5/serial:${CPLUS_INCLUDE_PATH}"
  export LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:${LIBRARY_PATH}"
  export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH}"
fi
