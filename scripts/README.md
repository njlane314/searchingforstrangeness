# Scripts

## Run 1 NuMI FHC original SAM definitions (pre-goodruns)

These are the original Run 1 NuMI FHC SAM definitions (pre-goodruns) you can
use as the persistent source when carrying filenames through processing
stages. This list includes beam, detector variations, strangeness (and
strangeness detvars), dirt, and EXT/data.

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

### Strangeness + detvars
- `prod_strange_resample_fhc_run1_fhc_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2`
- `Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2`
- `Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_wiremodX_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_wiremodYZ_sce_reco2_reco2`
- `Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2`
- `Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2`
- `detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2`

### Dirt
- `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0`
- `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1`

### EXT & Data (Run 1 FHC)
- `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2`
- `prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2` (already good; no extra
  good-runs pass needed)

## Next steps

1. Use the list above as the source definitions when you build new SAM defs.
   For example, to re-apply a good-runs condition (or a custom condition),
   run:

   ```bash
   ./scripts/apply_goodruns.sh <source_def> <goodruns_def> [condition]
   ```

2. To apply the good-runs condition to every Run 1 NuMI FHC definition above
   using consistent output names, run:

   ```bash
   ./scripts/apply_goodruns_run1_fhc.sh [--condition "<expr>"] [--dry-run]
   ```

3. If you need to split a nominal detector-variation sample into training
   and nominal subsets, use:

   ```bash
   ./scripts/split_detvar_stride.sh <source_def> <training_def> <nominal_def>
   ```

4. If you have multiple samples to split, create a batch file with whitespace
   triples and run:

   ```bash
   ./scripts/split_detvar_stride.sh --batch detvar_triples.txt
   ```

## End-to-end workflow: partition large samples, apply good-runs, and create
## training/validation splits

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

# Check the chunk sizes
samweb list-files --summary "defname: nl_ext_run1_chunk_0000_2500"
samweb list-files --summary "defname: nl_ext_run1_chunk_2500_2500"
samweb list-files --summary "defname: nl_ext_run1_chunk_5000_2500"
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
```

Name format used above:
- `<source_def>_goodruns` (or `_goodruns_reco2` when the source ends in
  `_good_reco2`)

### 3) Create training/validation splits from good-runs definitions

Once the good-runs definitions exist, create training/validation splits by
selecting `limit`/`offset` ranges. Example:

```bash
samweb create-definition nl_ext_run1_chunk_0000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with limit 2000"

samweb create-definition nl_ext_run1_chunk_0000_2500_val_500 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with offset 2000 with limit 500"
```

Name format used above:
- `<source_def>_train_<size>`
- `<source_def>_val_<size>`

Repeat steps 1–3 for other large sources (e.g., the beam-good Run 1 sample).

### Concrete, consistent naming commands for the Run 1 EXT and beam-good files

The commands below follow concise naming patterns and keep names consistent
across chunking, good-runs, and training/validation splits.

#### EXT (beam-off) sample

```bash
# Chunk the EXT sample (adjust size/offset as needed)
samweb create-definition nl_ext_run1_chunk_0000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with limit 2500"
samweb create-definition nl_ext_run1_chunk_2500_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 2500 with limit 2500"
samweb create-definition nl_ext_run1_chunk_5000_2500 \
  "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with offset 5000 with limit 2500"

# Apply good-runs to each chunk
./scripts/apply_goodruns.sh nl_ext_run1_chunk_0000_2500 \
  nl_ext_run1_chunk_0000_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_2500_2500 \
  nl_ext_run1_chunk_2500_2500_goodruns
./scripts/apply_goodruns.sh nl_ext_run1_chunk_5000_2500 \
  nl_ext_run1_chunk_5000_2500_goodruns

# Training/validation splits from good-runs chunks
samweb create-definition nl_ext_run1_chunk_0000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_0000_2500_val_500 \
  "defname: nl_ext_run1_chunk_0000_2500_goodruns with offset 2000 with limit 500"
samweb create-definition nl_ext_run1_chunk_2500_2500_train_2000 \
  "defname: nl_ext_run1_chunk_2500_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_2500_2500_val_500 \
  "defname: nl_ext_run1_chunk_2500_2500_goodruns with offset 2000 with limit 500"
samweb create-definition nl_ext_run1_chunk_5000_2500_train_2000 \
  "defname: nl_ext_run1_chunk_5000_2500_goodruns with limit 2000"
samweb create-definition nl_ext_run1_chunk_5000_2500_val_500 \
  "defname: nl_ext_run1_chunk_5000_2500_goodruns with offset 2000 with limit 500"
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

# Training/validation splits from good-runs chunks
samweb create-definition nl_beam_run1_chunk_0000_2500_train_2000 \
  "defname: nl_beam_run1_chunk_0000_2500_goodruns with limit 2000"
samweb create-definition nl_beam_run1_chunk_0000_2500_val_500 \
  "defname: nl_beam_run1_chunk_0000_2500_goodruns with offset 2000 with limit 500"
```

## Good-runs commands for every Run 1 NuMI FHC source above

Run the commands below to create a good-runs definition for each source listed
in this README. Beam and dirt keep the `<source_def>_goodruns` pattern; beam
detvars and strangeness detvars use shorter, standardized output names. The
beam-good data sample is already filtered, so it is intentionally omitted.

```bash
# Beam
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0 \
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1 \
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2 \
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3 \
  prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3_goodruns

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

# Strangeness + detvars
./scripts/apply_goodruns.sh prod_strange_resample_fhc_run1_fhc_reco2_reco2 \
  nl_run1_fhc_strange_nominal_goodruns
./scripts/apply_goodruns.sh detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2 \
  nl_run1_fhc_strange_cv_goodruns
./scripts/apply_goodruns.sh Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2 \
  nl_run1_fhc_strange_ly_down_goodruns
./scripts/apply_goodruns.sh Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2 \
  nl_run1_fhc_strange_ly_rayleigh_goodruns
./scripts/apply_goodruns.sh detvar_prod_strange_resample_fhc_run1_respin_wiremodX_sce_reco2_reco2 \
  nl_run1_fhc_strange_wiremodx_sce_goodruns
./scripts/apply_goodruns.sh detvar_prod_strange_resample_fhc_run1_respin_wiremodYZ_sce_reco2_reco2 \
  nl_run1_fhc_strange_wiremodyz_sce_goodruns
./scripts/apply_goodruns.sh Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2 \
  nl_run1_fhc_strange_wiremodyz_goodruns
./scripts/apply_goodruns.sh Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2 \
  nl_run1_fhc_strange_wiremodthetayz_goodruns
./scripts/apply_goodruns.sh detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2 \
  nl_run1_fhc_strange_sce_goodruns
./scripts/apply_goodruns.sh detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2 \
  nl_run1_fhc_strange_recomb2_goodruns

# Dirt
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0 \
  prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0_goodruns
./scripts/apply_goodruns.sh prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1 \
  prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1_goodruns

# EXT (Run 1 FHC)
./scripts/apply_goodruns.sh prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 \
  prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2_goodruns
```

## Expected outputs (good-runs + training/validation)

### Chunk definitions created in the examples above

```
# EXT chunks
nl_ext_run1_chunk_0000_2500
nl_ext_run1_chunk_2500_2500
nl_ext_run1_chunk_5000_2500

# Beam-good chunks
nl_beam_run1_chunk_0000_2500
nl_beam_run1_chunk_2500_2500
```

### Good-runs definitions created by the commands above

```
# Beam
prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0_goodruns
prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1_goodruns
prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2_goodruns
prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3_goodruns

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

# Strangeness + detvars
nl_run1_fhc_strange_nominal_goodruns
nl_run1_fhc_strange_cv_goodruns
nl_run1_fhc_strange_ly_down_goodruns
nl_run1_fhc_strange_ly_rayleigh_goodruns
nl_run1_fhc_strange_wiremodx_sce_goodruns
nl_run1_fhc_strange_wiremodyz_sce_goodruns
nl_run1_fhc_strange_wiremodyz_goodruns
nl_run1_fhc_strange_wiremodthetayz_goodruns
nl_run1_fhc_strange_sce_goodruns
nl_run1_fhc_strange_recomb2_goodruns

# Dirt
prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0_goodruns
prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1_goodruns

# EXT (Run 1 FHC)
prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2_goodruns
```

### Training/validation samples from the concise chunk workflow

When you follow the chunk + good-runs examples above (using 2500/2000/500), the
training and validation definitions created are:

```
# EXT chunks
nl_ext_run1_chunk_0000_2500_train_2000
nl_ext_run1_chunk_0000_2500_val_500
nl_ext_run1_chunk_2500_2500_train_2000
nl_ext_run1_chunk_2500_2500_val_500
nl_ext_run1_chunk_5000_2500_train_2000
nl_ext_run1_chunk_5000_2500_val_500

# Beam-good chunks
nl_beam_run1_chunk_0000_2500_train_2000
nl_beam_run1_chunk_0000_2500_val_500
```

If you create additional chunks, apply the same suffix pattern:
`<chunk_def>_train_<size>` and `<chunk_def>_val_<size>`.
