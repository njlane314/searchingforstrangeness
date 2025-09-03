# Script Guide

## Testing & Validation
- **check_containers.sh** – start an Apptainer container and verify Python, PyTorch, and MinkowskiEngine versions
- **eventdump.sh** – pull one file from a SAM definition and dump a single event with `lar`
- **inspect_root.sh** – list TTrees or top-level keys in a ROOT file using `uproot`
- **process.sh** – run a small sample through a FHiCL job and merge outputs for quick checks

## Job Workflow & Management
- **condor_lar.sh** – wrapper for submitting `lar` jobs with many configuration options
- **submit.sh** – iterate over workflow stages, calling `project.py` to clean, submit, and optionally check each stage
- **prestage.sh** – prestage predefined “strange” or “beam” datasets via `samweb`
- **summary.sh** – display SAM dataset summaries for default or user-specified definitions
- **project.py** – Python driver for project configuration, status tracking, and submission logic
- **jobstatus.sh** – report job IDs from logs or list current user jobs with monitoring links
- **render_workflow.py** – render workflow XML by applying entity mappings to stage templates

## Environment Setup & Packaging
- **init_payload_paths.sh** – export environment variables for FHiCL, weights, models, scripts, and adjust paths
- **pack_assets.sh** – bundle weights, models, scripts, and Python directories into a tarball, optionally upload to dropbox
- **tar_uboone.sh** – create a tarball from an MRB install or working directory after validating inputs
- **tar.sh** – call `tar_uboone.sh` and `pack_assets.sh` to produce repository and asset tarballs, exporting their paths
