#!/bin/sh
# run_fhicl_analysis.sh

# USAGE: source run_fhicl_analysis.sh
# Ensure you update the necessary input fields below

# ---------------------------------------------------------------------------------

# input: update these variables as per your case

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: source run_fhicl_analysis.sh <fhiclfile> <num_files>"
    exit 1
fi

# Assign input arguments to variables
fhiclfile=$1
num_files=$2

fhicl_base=$(basename "$fhiclfile" .fcl | sed 's/^run_//')

# Name of the SAM definition to query
samdef=prod_strange_resample_fhc_run2_fhc_reco2_reco2
#fhiclfile=/exp/uboone/app/users/nlane/production/KaonShortProduction01/srcs/ubana/ubana/searchingforstrangeness/run_emptyselectionfilter.fcl
output_directory="/exp/uboone/data/users/nlane/analysis"
combined_output="${output_directory}/${samdef}_${fhicl_base}_${num_files}_new_analysis.root"
tempdir="${output_directory}/temp_root_files"
#num_files=1

mkdir -p $tempdir

# ---------------------------------------------------------------------------------
BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
DEFAULT="\033[0m"

echo -e "${BLUE}Starting process for SAM definition: $samdef${DEFAULT}"

# ---------------------------------------------------------------------------------
# Get the first files from the SAMweb definition
# ---------------------------------------------------------------------------------

echo -e "${BLUE}Fetching the first $num_files files from SAM definition: $samdef${DEFAULT}"

files=$(samweb list-files defname:$samdef | head -n $num_files)

if [ -z "$files" ]; then
    echo -e "${RED}No files found in the SAM definition! Exiting...${DEFAULT}"
    exit 1
fi

echo -e "${YELLOW}Files fetched:${DEFAULT}"
echo "$files"

# ---------------------------------------------------------------------------------
# Run the FHiCL file on each of the files
# ---------------------------------------------------------------------------------

echo -e "${BLUE}Running the FHiCL file on the selected files...${DEFAULT}"

counter=0
for file in $files; do
    echo -e "${BLUE}Processing file: $file${DEFAULT}"

    # Locate the directory using samweb
    filedir=$(samweb locate-file $file | grep -o '/pnfs/.*' | head -n 1)

    if [ -z "$filedir" ]; then
        echo -e "${RED}Could not locate directory for file: $file. Skipping...${DEFAULT}"
        continue
    fi

    # Append the file name to the directory path
    filepath="${filedir}/${file}"

    if [ ! -f "$filepath" ]; then
        echo -e "${RED}File not found at: $filepath. Skipping...${DEFAULT}"
        continue
    fi

    # Output file name for this iteration
    outputfile="$tempdir/output_${counter}.root"

    echo -e "${BLUE}Running: lar -c $fhiclfile -s $filepath -T $outputfile${DEFAULT}"
    lar -c $fhiclfile -s $filepath -T $outputfile

    if [ ! -f $outputfile ]; then
        echo -e "${RED}FHiCL processing failed for file: $file${DEFAULT}"
        continue
    fi

    echo -e "${YELLOW}FHiCL processing successful. Output: $outputfile${DEFAULT}"

    counter=$((counter + 1))
done

# ---------------------------------------------------------------------------------
# Combine the resulting ROOT files using hadd
# ---------------------------------------------------------------------------------

echo -e "${BLUE}Combining all the output ROOT files...${DEFAULT}"

outputfiles=$(ls $tempdir/*.root)

if [ -z "$outputfiles" ]; then
    echo -e "${RED}No output ROOT files found to combine! Exiting...${DEFAULT}"
    exit 1
fi

echo -e "${YELLOW}Combining files into: $combined_output${DEFAULT}"
hadd -f $combined_output $outputfiles

if [ -f $combined_output ]; then
    echo -e "${GREEN}Successfully combined ROOT files into: $combined_output${DEFAULT}"
else
    echo -e "${RED}Combining ROOT files failed!${DEFAULT}"
fi

# ---------------------------------------------------------------------------------
# Clean up temporary files if desired
# ---------------------------------------------------------------------------------
echo -e "${BLUE}Cleaning up temporary files...${DEFAULT}"
rm -r $tempdir

echo -e "${GREEN}Process complete!${DEFAULT}"