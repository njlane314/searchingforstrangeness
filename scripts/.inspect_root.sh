#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/file.root" >&2
    exit 1
fi

FILE_PATH="$1"

python3 << EOF
import uproot
import sys

file_path = "$FILE_PATH"

try:
    print(f"Inspecting file: {file_path}")
    with uproot.open(file_path) as f:
        tree_names = f.keys(filter_classname="TTree")
        if not tree_names:
            print("No TTrees found in this file.")
            print("Available top-level keys:", f.keys())
        else:
            print("--- TTrees Found ---")
            for name in tree_names:
                print(f"\n--- Analyzing TTree: {name} ---")
                f[name].show()

except FileNotFoundError:
    print(f"Error: File not found at {file_path}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
    sys.exit(1)
EOF