#!/bin/sh
# run_fhicl_analysis.sh

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: source run_fhicl_analysis.sh <fhiclfile> <num_files>"
    return 1
fi

fhiclfile=$1
num_files=$2

fhicl_base=$(basename "$fhiclfile" .fcl | sed 's/^run_//')

#samdef=nl_ext_numi_fhc_beamoff_run1_reco2_training_2500
#samdef=nl_numi_fhc_beam_run1_reco2_training_250

#samdef=nl_numi_fhc_ext_run1_reco2_6000
#samdef=nl_numi_fhc_beam_run1_reco2_1000

#samdef=nl_prodgeni_numi_uboone_overlay_rhc_training_1600
#samdef=nl_numi_fhc_beam_run1_reco2_validation_250
#samdef=nl_lambda_nohadrons_reco2_validation_2000
#samdef=make_k0signal_overlay_testing_nohadrons_reco2_reco2
samdef=nl_strange_numi_fhc_run2_reco2_validation_982
output_directory="/exp/uboone/data/users/$USER/analysis"
combined_output="${output_directory}/${samdef}_${fhicl_base}_${num_files}_new_analysis.root"
tempdir="${output_directory}/temp_root_files"

mkdir -p $tempdir

BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
DEFAULT="\033[0m"

echo -e "${BLUE}Starting process for SAM definition: $samdef${DEFAULT}"

echo -e "${BLUE}Fetching the first $num_files files from SAM definition: $samdef${DEFAULT}"

files=$(samweb list-files defname:$samdef | head -n $num_files)

if [ -z "$files" ]; then
    echo -e "${RED}No files found in the SAM definition! Exiting...${DEFAULT}"
    return 1
fi

echo -e "${YELLOW}Files fetched:${DEFAULT}"
echo "$files"

echo -e "${BLUE}Running the FHiCL file on the selected files...${DEFAULT}"

counter=0
for file in $files; do
    echo -e "${BLUE}Processing file: $file${DEFAULT}"

    filedir=$(samweb locate-file $file | grep -o '/pnfs/.*' | sed 's/(\([0-9]*@[a-z0-9]*\))//g' | head -n 1)

    if [ -z "$filedir" ]; then
        echo -e "${RED}Could not locate directory for file: $file. Skipping...${DEFAULT}"
        continue
    fi

    filepath="${filedir}/${file}"

    if [ ! -f "$filepath" ]; then
        echo -e "${RED}File not found at: $filepath. Skipping...${DEFAULT}"
        continue
    fi

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


echo -e "${BLUE}Combining all the output ROOT files...${DEFAULT}"

outputfiles=$(ls $tempdir/*.root)

if [ -z "$outputfiles" ]; then
    echo -e "${RED}No output ROOT files found to combine! Exiting...${DEFAULT}"
    return 1
fi

echo -e "${YELLOW}Combining files into: $combined_output${DEFAULT}"
hadd -f $combined_output $outputfiles

if [ -f $combined_output ]; then
    echo -e "${GREEN}Successfully combined ROOT files into: $combined_output${DEFAULT}"
else
    echo -e "${RED}Combining ROOT files failed!${DEFAULT}"
fi


echo -e "${BLUE}Cleaning up temporary files...${DEFAULT}"
rm -r $tempdir

echo -e "${GREEN}Process complete!${DEFAULT}"
