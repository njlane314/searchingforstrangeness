#!/bin/bash

# Merge csv files from grid jobs together
# Useage: merge_csv.sh /path/to/training/data/folder/ 
# Execute in directory below the where all of your csv files are
# Training data folder needs to be somewhere on /exp/uboone/data/

destination=${1}

for plane in U V W; do

    echo "Plane $plane"
 
    csv_list="$(find ~+ -type f -name *_${plane}.csv)"

    rm ${destination}/training_output_${plane}.csv
    touch ${destination}/training_output_${plane}.csv

    while IFS= read -r line;
        do 
        echo "Adding $line"
        cat $line >> ${destination}/training_output_${plane}.csv
    done <<< "$csv_list"

done 
