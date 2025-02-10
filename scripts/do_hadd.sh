#!/bin/bash

outdir=/exp/uboone/data/users/nlane/strangeness/ana/
mkdir -p "${outdir}"

treename="${outdir}nlane_prod_strange_resample_fhc_run2_fhc_reco2_reco2_trainingimage_background_lambdamuon_ana.root"

rootdir=/pnfs/uboone/scratch/users/nlane/kaon_dl/v08_00_00_82/nlane_prod_strange_resample_fhc_run2_fhc_reco2_reco2_trainingimage_background_lambdamuon/out/

onelinefilelist=$(find "$rootdir" -type f -name "*.root" | tr '\n' ' ')

echo
echo "Found .root files:"
echo "$onelinefilelist"
echo

hadd -f "${treename}" ${onelinefilelist}
