#!/bin/bash
INPUT_DIR=$HOME/xenbase-gaf-pipeline/input-files

GOA_DIR=$INPUT_DIR/goa-gafs
NCBI_DIR=$INPUT_DIR/ncbi-maps
XB_DIR=$INPUT_DIR/xenbase-files

for DIR in "$GOA_DIR" "$NCBI_DIR" "$XB_DIR"; do
    echo "Scanning directory: $DIR"
    for FILENAME in "$DIR"/*; do
        if [ -f "$FILENAME" ]; then
            # Check if file is valid UTF-8
            if ! iconv -f UTF-8 -t UTF-8 "$FILENAME" -o /dev/null 2>/dev/null; then
                echo "$FILENAME -> File is NOT valid UTF-8!" 
            fi
        fi
    done
done

