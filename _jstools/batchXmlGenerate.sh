#!/bin/bash

if [ ! -f "../@jobs/cfgs/example.cfg" ]; then
    exit 1
fi

while read -r proj inputdef fcl numjobs; do
    [[ -z "$proj" || "$proj" =~ ^# ]] && continue
    ./xmlGenerate.sh -r v10_04_07_04 -t mc -u "$USER" -n "$proj" -i "$inputdef" -f "$fcl" -j "$numjobs" -o "../@jobs/$proj.xml"
done < projects.txt