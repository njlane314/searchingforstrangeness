#!/usr/bin/env bash

infile=${1:-definitions.txt}

if [[ ! -f "$infile" ]]; then
  echo "Input file '$infile' not found."
  exit 1
fi

while IFS= read -r def; do
  [[ -z "$def" ]] && continue
  [[ "$def" =~ ^# ]] && continue

  echo "================================================================"
  echo "Definition: $def"

  desc=$(samweb describe-definition "$def" 2>/dev/null)

  if [[ -z "$desc" ]]; then
    echo "  WARNING: samweb describe-definition failed or returned nothing."
  else
    created=$(echo "$desc" | sed -n 's/^ *Create date: *//p')
    owner=$(  echo "$desc" | sed -n 's/^ *Username: *//p')
    group=$(  echo "$desc" | sed -n 's/^ *Group: *//p')
    nfiles=$( echo "$desc" | sed -n 's/^ *Files: *//p')
    dims=$(   echo "$desc" | sed -n 's/^ *Dimensions: *//p')

    stage=$( echo "$dims" | sed -n "s/.*ub_project.stage *'\\([^']*\\)'.*/\\1/p")
    [[ -z "$stage" ]] && stage=$(echo "$dims" | sed -n "s/.*ub_project.stage *\\([^ ]*\\).*/\\1/p")

    version=$( echo "$dims" | sed -n "s/.*ub_project.version *'\\([^']*\\)'.*/\\1/p")
    [[ -z "$version" ]] && version=$(echo "$dims" | sed -n "s/.*ub_project.version *\\([^ ]*\\).*/\\1/p")

    projname=$( echo "$dims" | sed -n "s/.*ub_project.name *'\\([^']*\\)'.*/\\1/p")

    echo "  Definition summary:"
    [[ -n "$created" ]] && echo "    Created : $created"
    [[ -n "$owner"   ]] && echo "    Owner   : $owner"
    [[ -n "$group"   ]] && echo "    Group   : $group"
    [[ -n "$nfiles"  ]] && echo "    Files   : $nfiles"

    [[ -n "$projname" ]] && echo "    ub_project.name    : $projname"
    [[ -n "$stage"    ]] && echo "    ub_project.stage   : $stage"
    [[ -n "$version"  ]] && echo "    ub_project.version : $version"

    echo "    Dimensions:"
    echo "      $dims"
  fi

  echo
  echo "  Representative file meta-data:"

  rep_file=$(samweb list-definition-files "$def" 2>/dev/null | head -n 1)

  if [[ -z "$rep_file" ]]; then
    echo "    No files found (or samweb list-definition-files failed)."
    echo
    continue
  fi

  echo "    File: $rep_file"

  fmeta=$(samweb get-metadata "$rep_file" 2>/dev/null)

  if [[ -z "$fmeta" ]]; then
    echo "    WARNING: samweb get-metadata failed for this file."
    echo
    continue
  fi

  data_tier=$(echo "$fmeta" | sed -n 's/^ *data_tier *: *//p')
  file_type=$(echo "$fmeta" | sed -n 's/^ *file_type *: *//p')

  proj_name=$(        echo "$fmeta" | sed -n 's/^ *ub_project.name *: *//p')
  proj_stage=$(       echo "$fmeta" | sed -n 's/^ *ub_project.stage *: *//p')
  proj_version=$(     echo "$fmeta" | sed -n 's/^ *ub_project.version *: *//p')
  proj_release_ver=$( echo "$fmeta" | sed -n 's/^ *ub_project.release_version *: *//p')
  proj_production=$(  echo "$fmeta" | sed -n 's/^ *ub_project.production *: *//p')
  proj_campaign=$(    echo "$fmeta" | sed -n 's/^ *ub_project.campaign *: *//p')
  proj_release_type=$(echo "$fmeta" | sed -n 's/^ *ub_project.release_type *: *//p')

  app_family=$(  echo "$fmeta" | sed -n 's/^ *application.family *: *//p')
  app_name=$(    echo "$fmeta" | sed -n 's/^ *application.name *: *//p')
  app_version=$( echo "$fmeta" | sed -n 's/^ *application.version *: *//p')

  [[ -n "$data_tier" ]]        && echo "    data_tier              : $data_tier"
  [[ -n "$file_type" ]]        && echo "    file_type              : $file_type"

  [[ -n "$proj_name" ]]        && echo "    ub_project.name        : $proj_name"
  [[ -n "$proj_stage" ]]       && echo "    ub_project.stage       : $proj_stage"
  [[ -n "$proj_version" ]]     && echo "    ub_project.version     : $proj_version"
  [[ -n "$proj_release_ver" ]] && echo "    ub_project.release_version : $proj_release_ver"
  [[ -n "$proj_production" ]]  && echo "    ub_project.production  : $proj_production"
  [[ -n "$proj_campaign" ]]    && echo "    ub_project.campaign    : $proj_campaign"
  [[ -n "$proj_release_type" ]]&& echo "    ub_project.release_type: $proj_release_type"

  [[ -n "$app_family" ]]       && echo "    application.family     : $app_family"
  [[ -n "$app_name" ]]         && echo "    application.name       : $app_name"
  [[ -n "$app_version" ]]      && echo "    application.version    : $app_version"

  echo

done < "$infile"
