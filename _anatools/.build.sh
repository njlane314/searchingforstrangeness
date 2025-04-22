#!/bin/bash
#chmod +x build.sh

BIN_DIR="./bin"
mkdir -p "$BIN_DIR"

MACROS=(
    "plot_neutrino_energy.c plot_neutrino_energy"
    "plot_random_display.c plot_random_display"
    "plot_slice_info.c plot_slice_info"
    "plot_clarity_metrics.c plot_clarity_metrics"
)

COMPILE_CMD="g++ -o"

for macro in "${MACROS[@]}"; do
    read -r SOURCE EXEC <<< "$macro"
    $COMPILE_CMD "$EXEC" "$SOURCE" $(root-config --cflags --libs)
    if [ $? -eq 0 ]; then
        mv "$EXEC" "$BIN_DIR/"
    fi
done