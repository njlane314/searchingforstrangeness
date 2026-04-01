# Scripts

## Reco2 campaign source SAM definitions

These are the official Reco2 source definitions used by the source-based
campaign XML templates:

- `xml/numi_reco2_run1_fhc_campaign.xml`
- `xml/numi_reco2_run2a_fhc_campaign.xml`
- `xml/numi_reco2_run2b_rhc_campaign.xml`
- `xml/numi_reco2_run3b_rhc_campaign.xml`

Only the campaign-relevant categories are listed here: beam, dirt,
strangeness, EXT, and beam-on data. The full detector-variation inventories
remain in the internal note.

### Run 1 FHC
- Beam: `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0`,
  `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1`,
  `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2`,
  `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3`
- Dirt: `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0`,
  `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1`
- Strangeness: `prod_strange_resample_fhc_run1_fhc_reco2_reco2`
- EXT: `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2`
- Data: `prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2`

### Run 2a FHC
- Beam: `prodgenie_numi_overlay_detvar_CV_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2`
- Dirt: `prodgenie_numi_fhc_dirt_overlay_pandora_reco2_run2_reco2`
- Strangeness: `prod_strange_resample_fhc_run2_fhc_reco2_reco2`
- EXT: `prod_extnumi_swizzle_inclusive_v4_run2a_reco2_run2a_all_reco2`
- Data: `prod_numi_swizzle_inclusive_v4_run2_reco2_run2a_beam_good_reco2`

### Run 2b RHC
- Beam: `run2_numi_nu_overlay_pandora_unified_reco2_run2b_rhc_reco2`
- Dirt: `prodgenie_numi_rhc_dirt_overlay_pandora_reco2_run2_reco2`
- Strangeness: `prod_strange_resample_rhc_run2_rhc_reco2_reco2`
- EXT: `prod_extnumi_swizzle_crt_inclusive_v4b_offbeam_run2_reco2_new2_run2_reco2_all`
- Data: `prod_numi_swizzle_inclusive_v4b_run2_beamon_run2_reco2_run2_beam_good_reco2`

### Run 3b RHC
- Beam: `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0`,
  `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample1`,
  `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample2_v3`,
  `prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_sample3`
- Dirt: `prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_sample0`,
  `prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_sample1`
- Strangeness: `prod_strange_resample_rhc_run3_rhc_reco2_reco2`
- EXT: `prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2`
- Data: `prod_numi_mcc9_v08_00_00_45_run3b_run3b_reco2_beam_good_reco2`

Run 1 RHC is not templated here because this repository does not record a
matching dedicated Run 1 RHC strangeness source definition.

## Run 1 NuMI FHC source SAM definitions

These are the Run 1 NuMI FHC SAM definitions to use as the persistent source
when carrying filenames through processing stages. This list includes beam,
detector variations, strangeness (and strangeness detvars), dirt, and
EXT/data. The checked-in XML workflows use these source definitions directly.

### Beam
- `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0`
- `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1`
- `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2`
- `prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3`

### Detector variations (Run 1 FHC)
- `prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2`
- `prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2`

### Strangeness (nominal)
- `prod_strange_resample_fhc_run1_fhc_reco2_reco2`

### Strangeness detvars confirmed in the Run 1 FHC Reco2 campaign tracking
- `detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2`
- `Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2`
- `Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run_respin_wiremodX_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run_respin_wiremodYZ_sce_reco2_reco2`
- `Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2`
- `Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2`

Notes:
- No separate Run 1 FHC strangeness Reco2 LY-attenuation sample is currently
  listed in the confirmed production table.
- One production row labels
  `Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2`
  as "WireMod ThetaXZ", but the SAM definition name itself says `WireMod_YZ`.
  This repository keeps the literal SAM definition until that campaign label is
  resolved externally.

### Dirt
- `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0`
- `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1`

### EXT & Data (Run 1 FHC)
- `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2`
- `prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2` (already good; no extra
  good-runs pass needed)

## Next steps

1. Use the list above as the source definitions when you build new SAM defs.
   If you still need the legacy good-runs filtering helpers for an ad hoc SAM
   workflow, run:

   ```bash
   ./scripts/apply_goodruns.sh <source_def> <goodruns_def> [condition]
   ```

2. To apply the legacy good-runs condition to every Run 1 NuMI FHC definition
   above using consistent output names, run:

   ```bash
   ./scripts/apply_goodruns_run1_fhc.sh [--condition "<expr>"] [--dry-run]
   ```

3. If you need to split a nominal detector-variation sample into orthogonal
   training and plotting (nominal) subsets, use:

   ```bash
   ./scripts/split_detvar_stride.sh <source_def> <training_def> <plotting_def>
   ```

   The stride-2 split guarantees that training and plotting samples are
   orthogonal.

4. If you have multiple samples to split, create a batch file with whitespace
   triples and run:

   ```bash
   ./scripts/split_detvar_stride.sh --batch detvar_triples.txt
   ```

5. To process training definitions on the grid without running inference, use
   `xml/numi_fhc_run1_training.xml` and update its input SAM definitions and
   `numjobs` counts to match the training definitions you created.

6. To process the stride-split nominal detvar training definitions created by
   `split_detvar_stride.sh`, use `xml/numi_fhc_run1_detvar_training.xml`.
   That workflow is image-plus-selection only and uses the no-inference
   selection wrapper.

## Legacy end-to-end workflow: partition large samples, apply good-runs, and
## create training splits

This workflow is intended for Run 1 NuMI FHC sources in the list above. It
gives consistent naming, and keeps the steps repeatable.

### 1) Partition large sources into smaller chunks

Pick a chunk size (e.g., 2500 files). Use `limit`/`offset` to define chunks.
Example for the EXT sample (repeat for any large source):

```bash
# Chunk 1
samweb create-definition nl_ext_run1_chunk_0000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with limit 2500"

# Chunk 2
samweb create-definition nl_ext_run1_chunk_2500_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 2500 with limit 2500"

# Chunk 3
samweb create-definition nl_ext_run1_chunk_5000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 5000 with limit 2500"

# Chunk 4
samweb create-definition nl_ext_run1_chunk_7500_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 7500 with limit 2500"

# Check the chunk sizes
samweb list-files --summary "defname: nl_ext_run1_chunk_0000_2500"
samweb list-files --summary "defname: nl_ext_run1_chunk_2500_2500"
samweb list-files --summary "defname: nl_ext_run1_chunk_5000_2500"
samweb list-files --summary "defname: nl_ext_run1_chunk_7500_2500"
```

Name format used above:
- `nl_<tag>_chunk_<offset>_<size>`

### 2) Apply good-runs with consistent output names

Apply good-runs to each chunk (or to the full definitions if you skip chunking):

```bash
./scripts/apply_goodruns.sh nl_ext_run1_chunk_0000_2500 \
  nl_ext_run1_chunk_0000_2500_goodruns

./scripts/apply_goodruns.sh nl_ext_run1_chunk_2500_2500 \
  nl_ext_run1_chunk_2500_2500_goodruns

./scripts/apply_goodruns.sh nl_ext_run1_chunk_5000_2500 \
  nl_ext_run1_chunk_5000_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_7500_2500 \
  nl_ext_run1_chunk_7500_2500_goodruns
```

Name format used above:
- `<source_def>_goodruns` (or `_goodruns_reco2` when the source ends in
  `_good_reco2`)

### 3) Create training splits from good-runs definitions

Once the good-runs definitions exist, create training splits that are disjoint
from any plotting samples by selecting `limit`/`offset` ranges. Example:

```bash
samweb create-definition nl_ext_run1_chunk_0000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with limit 2000"

```

Name format used above:
- `<source_def>_train_<size>`

Repeat steps 1–3 for other large sources (e.g., the beam-good Run 1 sample).

### Concrete, consistent naming commands for the Run 1 EXT and beam-good files

The commands below follow concise naming patterns and keep names consistent
across chunking, good-runs, and training splits.

#### EXT (beam-off) sample

```bash
# Chunk the EXT sample (adjust size/offset as needed)
samweb create-definition nl_ext_run1_chunk_0000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with limit 2500"
samweb create-definition nl_ext_run1_chunk_2500_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 2500 with limit 2500"
samweb create-definition nl_ext_run1_chunk_5000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 5000 with limit 2500"
samweb create-definition nl_ext_run1_chunk_7500_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 7500 with limit 2500"

# Apply good-runs to each chunk
./scripts/apply_goodruns.sh nl_ext_run1_chunk_0000_2500 \
  nl_ext_run1_chunk_0000_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_2500_2500 \
  nl_ext_run1_chunk_2500_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_5000_2500 \
  nl_ext_run1_chunk_5000_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_7500_2500 \
  nl_ext_run1_chunk_7500_2500_goodruns

# Training splits from good-runs chunks
samweb create-definition nl_ext_run1_chunk_0000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_2500_2500_train_2000 \
  "defname: nl_ext_run1_chunk_2500_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_5000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_5000_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_7500_2500_train_2000 \
  "defname: nl_ext_run1_chunk_7500_2500_goodruns with limit 2000"
```

#### Beam-good sample

```bash
# Chunk the beam-good sample (adjust size/offset as needed)
samweb create-definition nl_beam_run1_chunk_0000_2500 \
  "defname: prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2 with limit 2500"
samweb create-definition nl_beam_run1_chunk_2500_2500 \
  "defname: prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2 with offset 2500 with limit 2500"

# Apply good-runs (note the source already includes _good_reco2)
./scripts/apply_goodruns.sh nl_beam_run1_chunk_0000_2500 \
  nl_beam_run1_chunk_0000_2500_goodruns
./scripts/apply_goodruns.sh nl_beam_run1_chunk_2500_2500 \
  nl_beam_run1_chunk_2500_2500_goodruns

# Training splits from good-runs chunks
samweb create-definition nl_beam_run1_chunk_0000_2500_train_2000 \
  "defname: nl_beam_run1_chunk_0000_2500_goodruns with limit 2000"
```

## Legacy good-runs commands for every Run 1 NuMI FHC source above

Run the commands below to create a good-runs definition for the generic Run 1
FHC beam, detvar, dirt, and EXT sources. Beam, dirt, and EXT use the same
short `nl_run1_fhc_*` naming style as the detector variations for consistency.
The beam-good data sample is already filtered, so it is intentionally omitted.
The dedicated strangeness Reco2 SAM definitions are also intentionally omitted:
their upstream Reco1 source names already indicate that the good-runs
requirement was applied before the strangeness productions were made.

```bash
# Beam
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0 \
  nl_run1_fhc_beam_sample0_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1 \
  nl_run1_fhc_beam_sample1_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2 \
  nl_run1_fhc_beam_sample2_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3 \
  nl_run1_fhc_beam_sample3_goodruns

# Detector variations (Run 1 FHC)
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2 \
  nl_run1_fhc_detvar_cv_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_ly_supp75_att8m_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_ly_rayleigh_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_ly_down_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2 \
  nl_run1_fhc_detvar_sce_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_recomb2_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2 \
  nl_run1_fhc_detvar_wiremodx_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_wiremodyz_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2 \
  nl_run1_fhc_detvar_wiremodthetaxz_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2 \
  nl_run1_fhc_detvar_wiremodthetayz_goodruns

# Dirt
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0 \
  nl_run1_fhc_dirt_sample0_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1 \
  nl_run1_fhc_dirt_sample1_goodruns

# EXT (Run 1 FHC)
./scripts/apply_goodruns.sh prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 \
  nl_run1_fhc_ext_goodruns
```

## Expected outputs (good-runs + training)

### Chunk definitions created in the examples above

```
# EXT chunks
nl_ext_run1_chunk_0000_2500
nl_ext_run1_chunk_2500_2500
nl_ext_run1_chunk_5000_2500
nl_ext_run1_chunk_7500_2500

# Beam-good chunks
nl_beam_run1_chunk_0000_2500
nl_beam_run1_chunk_2500_2500
```

### Good-runs definitions created by the commands above

```
# Beam
nl_run1_fhc_beam_sample0_goodruns
nl_run1_fhc_beam_sample1_goodruns
nl_run1_fhc_beam_sample2_goodruns
nl_run1_fhc_beam_sample3_goodruns

# Detector variations (Run 1 FHC)
nl_run1_fhc_detvar_cv_goodruns
nl_run1_fhc_detvar_ly_supp75_att8m_goodruns
nl_run1_fhc_detvar_ly_rayleigh_goodruns
nl_run1_fhc_detvar_ly_down_goodruns
nl_run1_fhc_detvar_sce_goodruns
nl_run1_fhc_detvar_recomb2_goodruns
nl_run1_fhc_detvar_wiremodx_goodruns
nl_run1_fhc_detvar_wiremodyz_goodruns
nl_run1_fhc_detvar_wiremodthetaxz_goodruns
nl_run1_fhc_detvar_wiremodthetayz_goodruns

# Dirt
nl_run1_fhc_dirt_sample0_goodruns
nl_run1_fhc_dirt_sample1_goodruns

# EXT (Run 1 FHC)
nl_run1_fhc_ext_goodruns
```

### Training samples from the concise chunk workflow

When you follow the chunk + good-runs examples above (using 2500/2000), the
training definitions created are:

```
# EXT chunks
nl_ext_run1_chunk_0000_2500_train_2000
nl_ext_run1_chunk_2500_2500_train_2000
nl_ext_run1_chunk_5000_2500_train_2000
nl_ext_run1_chunk_7500_2500_train_2000

# Beam-good chunks
nl_beam_run1_chunk_0000_2500_train_2000
```

If you create additional chunks, apply the same suffix pattern:
`<chunk_def>_train_<size>`.
