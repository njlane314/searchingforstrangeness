# Scripts

## Current XML surface

`xml/` now contains these staged campaign entry points:

- `xml/numi_run1_fhc_campaign.xml`
- `xml/numi_run2a_fhc_campaign.xml`
- `xml/numi_run2b_rhc_campaign.xml`
- `xml/numi_run3b_rhc_campaign.xml`
- `xml/numi_run4a_rhc_campaign.xml`
- `xml/numi_run4b_rhc_campaign.xml`
- `xml/numi_run4c_fhc_campaign.xml`
- `xml/numi_run4d_fhc_campaign.xml`
- `xml/numi_run5_fhc_campaign.xml`

Legacy XMLs were moved to `reference/legacy_xml/` so they stay available as
reference without looking like active submission entry points.

The checked-in campaign XMLs use one stage per sample:

- G4-updated MC beam chains run `redk2nu -> evtw -> fullchain`
- legacy MC dirt chains still run `redk2nu -> evtw -> fullchain_oldflux_rw`
- dedicated strangeness MC chains already have the new NuMI flux, so they run
  `redk2nu -> evtw -> fullchain` without old-flux reweighting
- data and EXT chains run `fullchain_data`
- detector-variation chains run a single fullchain stage
- only `fcl_evtw_00` is active, so the MC chains keep the first 100
  multisim universes
- active campaign XMLs default to 25-job test submissions

For Run 1 specifically:

- `xml/numi_run1_fhc_campaign.xml` is the active Run 1 test campaign and only
  carries `beam`, `dirt`, `strangeness`, and `ext`
- the checked-in Run 1 EXT surface is still being treated as FHC-era by
  temporary convention
- that Run 1 test campaign skips inference and stops at
  `imageprod -> run_nuselection_*_slim`

The checked-in local validation helper now uses the staged EventWeight dev
wrappers directly, and those wrappers run redk2nu internally before the
standard EventWeight configuration.

## Source SAM definitions

These are the source definitions currently encoded in the staged campaign XMLs.

### Run 1 FHC

- Beam shards:
  - `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0`
  - `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1`
  - `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2`
  - `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3`
- Dirt shards:
  - `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0`
  - `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1`
- Dedicated strangeness:
  - `prod_strange_resample_fhc_run1_fhc_reco2_reco2`
- EXT:
  - `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2`

### Run 2a FHC

- Beam:
  - `prodgenie_numi_overlay_detvar_CV_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
- Dirt:
  - `prodgenie_numi_fhc_dirt_overlay_pandora_reco2_run2_reco2`
- Dedicated strangeness:
  - `prod_strange_resample_fhc_run2_fhc_reco2_reco2`
- EXT:
  - `prod_extnumi_swizzle_inclusive_v4_run2a_reco2_run2a_all_reco2`
- Data:
  - `prod_numi_swizzle_inclusive_v4_run2_reco2_run2a_beam_good_reco2`
- Beam detvars:
  - `prodgenie_numi_nue_overlay_detvar_CV_run2_FHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_overlay_detvar_LYAttenuation_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_overlay_detvar_LYRayleigh_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_overlay_detvar_LYDown_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_SCE_run2_FHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_overlay_detvar_Recomb2_run2_FHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_detvar_WireModX_standard_FHC_nu_overlay_400k_run2_reco2_reco2_reco2`
  - `prodgenie_numi_overlay_detvar_WireModYZ_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_FHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_pandora_reco2_run2_FHC_WireModThetaYZ_withSplines_reco2`

### Run 2b RHC

- Beam:
  - `run2_numi_nu_overlay_pandora_unified_reco2_run2b_rhc_reco2`
- Dirt:
  - `prodgenie_numi_rhc_dirt_overlay_pandora_reco2_run2_reco2`
- Dedicated strangeness:
  - `prod_strange_resample_rhc_run2_rhc_reco2_reco2`
- EXT:
  - `prod_extnumi_swizzle_crt_inclusive_v4b_offbeam_run2_reco2_new2_run2_reco2_all`
- Data:
  - `prod_numi_swizzle_inclusive_v4b_run2_beamon_run2_reco2_run2_beam_good_reco2`
- Beam detvars:
  - `prodgenie_numi_nu_overlay_detvar_CV_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_LYAttenuation_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYRayleigh_run2_reco2`
  - `prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYDown_run2_reco2`
  - `prodgenie_numi_overlay_run2_RHC_standard_nu_reco2_detvar_SCE_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_Recomb2_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModX_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModYZ_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_RHC_reco2_v08_00_00_55_run2_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_wsplines_run2_RHC_reco2_v08_00_00_5_run2_reco2`

### Run 3b RHC

- Beam shards:
  - `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0`
  - `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample1`
  - `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample2_v3`
  - `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample3`
- Dirt shards:
  - `prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_sample0`
  - `prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_sample1`
- Dedicated strangeness:
  - `prod_strange_resample_rhc_run3_rhc_reco2_reco2`
- EXT:
  - `prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2`
- Data:
  - `prod_numi_mcc9_v08_00_00_45_run3b_run3b_reco2_beam_good_reco2`
- Beam detvars:
  - `numi_run3_standard_nu_overlay_cv_reco2_v08_00_00_54_run3_reco2`
  - `numi_run3_standard_nu_overlay_LYattenuation_reco2_run3_reco2`
  - `numi_run3_standard_nu_overlay_LYrayleigh_reco2_run3_reco2`
  - `numi_run3_standard_nu_overlay_LYDown_reco2_v08_00_00_54_run3_reco2`
  - `numi_run3_standard_nu_overlay_SCE_reco2_v08_00_00_54_run3_reco2`
  - `numi_run3_standard_nu_overlay_Recomb2_reco2_v08_00_00_54_run3_reco2`
  - `numi_run3_standard_nu_overlay_cv_v08_00_00_54_run3_reco2_WiremodX_run3b_reco2`
  - `numi_run3_standard_nu_overlay_WireModYZ_reco2_run3b_reco2`
  - `numi_run3_standard_nu_overlay_WireModThetaXZ_reco2_v08_00_00_54_run3_reco2`
  - `numi_run3_standard_nu_overlay_WireModThetaYZwithsplines_reco2_run3_reco2`

## Good-runs notes

The repo no longer keeps helper scripts for good-runs SAM-definition creation.
The explicit `samweb create-definition ...` commands now live at the bottom of
`SAMPLES`, alongside the `original reco2`, `goodruns list`, and `updated
reco2` mapping for each sample.


## Useful commands

Inspect the current source definitions and notes:

```bash
./describe-definitions.sh
```

Count the active campaign definitions and flag anything above the default
5000-file sharding threshold:

```bash
./scripts/campaign-jobs.sh
```

Run the checked-in dev FHiCL chain locally before you launch a campaign:

```bash
./scripts/validate-campaign.sh --workflow mc --samdef prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0
```

Run the Run 1 data/EXT-like staged path locally:

```bash
./scripts/validate-campaign.sh --workflow data
```

Build the compact local ntuple validation surface locally:

```bash
./scripts/validate-campaign.sh --workflow ntuple
```

Count a specific campaign XML:

```bash
./scripts/campaign-jobs.sh --xml xml/numi_run1_fhc_campaign.xml
```

Or use the wrapper that prints source counts and then applies the whole plan,
targeting roughly 100000 events in the first shard by default while keeping
the second shard untrimmed:

```bash
./scripts/train-template.sh
```

Preview the sharding commands without creating any SAM definitions:

```bash
./scripts/train-template.sh --dry-run
```

To keep the original full alternating-file shards through the wrapper:

```bash
./scripts/train-template.sh --full-shards
```
