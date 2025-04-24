#!/bin/bash
#chmod +x build.sh

BIN_DIR="./bin"
mkdir -p "$BIN_DIR"

MACROS=(
    "plot_random_display.c plot_random_display"
)

COMPILE_CMD="g++ -o"

for macro in "${MACROS[@]}"; do
    read -r SOURCE EXEC <<< "$macro"
    $COMPILE_CMD "$EXEC" "$SOURCE" $(root-config --cflags --libs)
    if [ $? -eq 0 ]; then
        mv "$EXEC" "$BIN_DIR/"
    fi
done