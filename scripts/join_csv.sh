#!/bin/bash

destination=${1}

if [[ -z "$destination" ]]; then
    echo "Error: No destination folder provided."
    exit 1
fi

temp_dir="./temp"
mkdir -p "$temp_dir"

echo "Extracting .csv files from log.tar files..."
find log -type f -name "log.tar" | while read tarfile; do
    prefix=$(basename "$(dirname "$tarfile")")
    echo "Processing $tarfile with prefix $prefix"

    tar -xf "$tarfile" --wildcards '*.csv'

    for file in ./*.csv; do
        if [[ -f "$file" ]]; then
            mv "$file" "$temp_dir/${prefix}_$(basename "$file")" && echo "Moved $file to $temp_dir/${prefix}_$(basename "$file")" || echo "Failed to move $file"
        fi
    done
done

for plane in U V W; do
    echo "Processing Plane $plane"

    csv_list=$(find "$temp_dir" -type f -name "*_${plane}.csv")

    if [[ -z "$csv_list" ]]; then
        echo "No files found for Plane $plane. Skipping."
        continue
    fi

    output_file="${destination}/training_output_${plane}.csv"
    rm -f "$output_file"
    touch "$output_file"

    for file in $csv_list; do
        echo "Adding $file to $output_file"
        cat "$file" >> "$output_file"
    done
done

echo "Merging completed. Merged files are in $destination."

