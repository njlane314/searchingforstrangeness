#!/bin/bash

outdir=/exp/uboone/data/users/nlane/analysis/
mkdir -p "${outdir}"

treename="${outdir}numi_fhc_overlay_intrinsic_strangeness_run1"

rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_strange_run1_eventselectionfilter/out/
#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_run1_eventselectionfilter/
#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_run1_eventselectionfilter/out/
#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_lambda_nohadrons_eventselectionfilter_trainingsample_10_test/out/
#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_lambda_nohadrons_eventselectionfilter_trainingsample_10_test/out/

onelinefilelist=$(find "$rootdir" -type f -name "*.root" | tr '\n' ' ')

echo
echo "Found .root files:"
echo "$onelinefilelist"
echo

hadd -f "${treename}" ${onelinefilelist}
