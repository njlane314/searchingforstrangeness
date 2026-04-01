# Scripts

## Current XML surface

`xml/` now contains only the four staged Reco2 campaign entry points:

- `xml/numi_reco2_run1_fhc_campaign.xml`
- `xml/numi_reco2_run2a_fhc_campaign.xml`
- `xml/numi_reco2_run2b_rhc_campaign.xml`
- `xml/numi_reco2_run3b_rhc_campaign.xml`

Legacy XMLs were moved to `reference/legacy_xml/` so they stay available as
reference without looking like active submission entry points.

All four checked-in campaign XMLs follow the same high-level rules:

- nominal beam, dirt, and dedicated strangeness chains run
  `redk2nu -> evtw -> image -> sel`
- only `fcl_evtw_00` is active, so the nominal MC chains keep the first 100
  multisim universes
- detector-variation chains run `image -> sel` only
- inference stages are omitted because the checked-in model bundle has not been
  updated yet
- `numjobs` values are placeholders and should be replaced with
  `samweb count-files defname:<...>` results before submission

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
- Data:
  - `prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2`
- Generic Run 1 beam detvars:
  - `prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2`
  - `prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2`

### Run 1 FHC dedicated strangeness detvars

These are the confirmed Run 1 FHC strangeness Reco2 detvars currently wired
into `xml/numi_reco2_run1_fhc_campaign.xml`:

- `Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2`
- `Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run_respin_wiremodX_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run_respin_wiremodYZ_sce_reco2_reco2`
- `Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2`
- `Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2`

Notes:

- No separate Run 1 FHC strangeness LY-attenuation Reco2 sample is currently
  confirmed here.
- One campaign-tracking row labels
  `Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2`
  as ThetaXZ, but the SAM definition itself says `WireMod_YZ`. The XML keeps
  that sample under a neutral `wiremod_mismatch` label until the bookkeeping
  is clarified externally.

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

Current repo policy is:

- generic overlay, dirt, generic detvars, and EXT can be replaced with
  upstream good-run-filtered definitions if that is the campaign surface you
  want to process
- dedicated strangeness Reco2 definitions should be used directly and should
  not be re-filtered here, because their upstream provenance already indicates
  good-run filtering
- beam-good data definitions should be used directly

The checked-in XMLs currently preserve the source definitions above; they do
not auto-rewrite anything to good-run-filtered variants.

## Run 1 CV train/template shards

Only Run 1 currently has the CV split encoded directly into the staged
campaign XML. The split happens at the SAM-definition level so training and
template surfaces are orthogonal before any processing stage runs.

The checked-in batch file is:

- `scripts/run1_detvar_cv_shards.txt`

Create the shard definitions with:

```bash
./scripts/split_detvar_stride.sh --batch scripts/run1_detvar_cv_shards.txt
```

That produces:

- `nl_run1_fhc_beam_detvar_cv_train_shard`
- `nl_run1_fhc_beam_detvar_cv_template_shard`
- `nl_run1_fhc_strangeness_detvar_cv_train_shard`
- `nl_run1_fhc_strangeness_detvar_cv_template_shard`

The Run 1 campaign XML already points at those names. If you use different SAM
definition names, update the four XML entities near the top of
`xml/numi_reco2_run1_fhc_campaign.xml`.

## Useful commands

Inspect the current source definitions and notes:

```bash
./describe-definitions.sh
```

Count the active campaign definitions and flag anything above the default
5000-file sharding threshold:

```bash
./scripts/query_campaign_numjobs.sh
```

Run the checked-in dev FHiCL chain locally before you launch a campaign:

```bash
./scripts/validate_campaign_local.sh --workflow mc --samdef prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0
```

Run the Run 1 data/EXT-like staged path locally:

```bash
./scripts/validate_campaign_local.sh --workflow data
```

Build the compact amarantin validation surface locally:

```bash
./scripts/validate_campaign_local.sh --workflow amarantin
```

Count a specific campaign XML and include derived train/template shard
definitions:

```bash
./scripts/query_campaign_numjobs.sh \
  --xml xml/numi_reco2_run1_fhc_campaign.xml \
  --include-derived-shards
```

Create one orthogonal pair by hand:

```bash
./scripts/split_detvar_stride.sh \
  prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2 \
  nl_run1_fhc_beam_detvar_cv_train_shard \
  nl_run1_fhc_beam_detvar_cv_template_shard
```

Create the whole Run 1 CV shard set:

```bash
./scripts/split_detvar_stride.sh --batch scripts/run1_detvar_cv_shards.txt
```

Or use the wrapper that prints source counts and then applies the whole plan:

```bash
./scripts/create_training_template_shards.sh
```

Preview the sharding commands without creating any SAM definitions:

```bash
./scripts/create_training_template_shards.sh --dry-run
```
