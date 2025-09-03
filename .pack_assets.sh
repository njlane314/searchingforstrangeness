#!/usr/bin/env bash
# pack_assets.sh — build an RCDS-ready tarball for grid jobs
# Usage:
#   ./pack_assets.sh -r <assets_root> -o <out.tar.gz> [-d <dropbox_dir>]
#
# Example:
#   ./pack_assets.sh \
#     -r /home/$USER/my-ana/assets \
#     -o /home/$USER/assets.tar.gz \
#     -d /dune/app/users/$USER
#
# Or for uboone:
#   ./pack_assets.sh -r /home/$USER/assets -o /home/$USER/assets.tar.gz -d /uboone/app/users/$USER
#
# Then submit with:
#   --tar_file_name=dropbox:///dune/app/users/$USER/assets.tar.gz
#   # or
#   --tar_file_name=dropbox:///uboone/app/users/$USER/assets.tar.gz

set -euo pipefail

usage() {
  echo "Usage: $0 -r <assets_root> -o <out.tar.gz> [-d <dropbox_dir>]" >&2
  echo "  <assets_root> should contain any of: weights/, models/, scripts/, python/" >&2
  exit 2
}

ASSETS_ROOT=""
OUT_TAR=""
DROPBOX_DIR=""

while getopts ":r:o:d:" opt; do
  case "$opt" in
    r) ASSETS_ROOT="$OPTARG" ;;
    o) OUT_TAR="$OPTARG" ;;
    d) DROPBOX_DIR="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "${ASSETS_ROOT}" || -z "${OUT_TAR}" ]] && usage
[[ -d "$ASSETS_ROOT" ]] || { echo "ERROR: assets root '$ASSETS_ROOT' not found"; exit 1; }

# Tools
command -v tar >/dev/null || { echo "ERROR: 'tar' not found"; exit 1; }
HAVE_PIGZ=0; command -v pigz >/dev/null && HAVE_PIGZ=1

# Make a tidy staging area with only what we need
STAGE="$(mktemp -d)"
cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

# Copy known subdirs if present; you can add more patterns here if needed
shopt -s nullglob dotglob
pushd "$ASSETS_ROOT" >/dev/null
for d in weights models scripts python; do
  if [[ -d "$d" ]]; then
    mkdir -p "$STAGE/$d"
    rsync -a --delete --chmod=F=u+rw,Da+rx \
      --exclude=".git" --exclude="__pycache__" --exclude="*.pyc" \
      --exclude=".DS_Store" --exclude=".ipynb_checkpoints" \
      "$d/" "$STAGE/$d/"
  fi
done

# Optional: include any top-level config files you want alongside the dirs
for f in *.json *.yaml *.yml *.txt; do
  [[ -e "$f" ]] && install -m 0644 "$f" "$STAGE/"
done
popd >/dev/null

# Create a simple manifest inside the tarball (file list + sha256)
if command -v sha256sum >/dev/null; then
  ( cd "$STAGE" && find . -type f -print0 | sort -z | xargs -0 sha256sum ) > "$STAGE/MANIFEST.sha256" || true
fi

# Reproducible tar flags (GNU tar): stable order, fixed owner, fixed mtime
TAR_FLAGS=(--numeric-owner --owner=0 --group=0 --sort=name --mtime='UTC 2025-01-01' -C "$STAGE" .)

# Compression
if (( HAVE_PIGZ )); then
  # fastest multi-core gzip
  tar -cf - "${TAR_FLAGS[@]}" | pigz -9 > "$OUT_TAR"
else
  tar -czf "$OUT_TAR" "${TAR_FLAGS[@]}"
fi

# Print size + checksum
echo "Wrote: $OUT_TAR"
du -h "$OUT_TAR" | awk '{print "Size: " $1}'
if command -v sha256sum >/dev/null; then
  sha256sum "$OUT_TAR" | awk '{print "SHA256: " $1}'
fi

# Optionally copy into RCDS dropbox
if [[ -n "$DROPBOX_DIR" ]]; then
  mkdir -p "$DROPBOX_DIR"
  cp -f "$OUT_TAR" "$DROPBOX_DIR/"
  BASENAME="$(basename "$OUT_TAR")"
  echo "Copied to: $DROPBOX_DIR/$BASENAME"
  echo
  echo "Use this in job submission:"
  echo "  --tar_file_name=dropbox://$DROPBOX_DIR/$BASENAME"
fi
