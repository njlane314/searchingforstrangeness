#!/usr/bin/env bash
# Print the largest directory and file in this repository.
# Usage: ./find_largest.sh
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel 2>/dev/null || pwd)

# Find largest directory
largest_dir_line=$(find "$repo_root" -mindepth 1 -type d -print0 | du --files0-from=- -sb 2>/dev/null | sort -nr | head -n 1)
dir_size=$(awk '{print $1}' <<<"$largest_dir_line")
dir_path=$(cut -f2 <<<"$largest_dir_line")
dir_human=$(numfmt --to=iec "$dir_size" 2>/dev/null)

# Find largest file
largest_file_line=$(find "$repo_root" -type f -printf '%s\t%p\n' | sort -nr | head -n 1)
file_size=$(awk '{print $1}' <<<"$largest_file_line")
file_path=$(cut -f2 <<<"$largest_file_line")
file_human=$(numfmt --to=iec "$file_size" 2>/dev/null)

echo "Largest directory: $dir_path ($dir_human)"
echo "Largest file: $file_path ($file_human)"

