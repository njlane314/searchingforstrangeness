#!/bin/bash

outdir=/exp/uboone/data/users/nlane/analysis/
mkdir -p "${outdir}"

treename="${outdir}numi_fhc_run1_strangeness_ana.root"

#treename="${outdir}numi_fhc_run1_detvar_recomb2_strangeness_ana.root"

#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_eventselectionfilter_detvar_recomb2_strangeness/out/
#rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_eventselectionfilter_numi_fhc_run1_beam/out/ 
rootdir=/pnfs/uboone/scratch/users/nlane/strangeness/v08_00_00_82/nl_numi_fhc_eventselectionfilter_numi_fhc_run1_strangeness/out/
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
