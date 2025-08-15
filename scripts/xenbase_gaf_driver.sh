#!/bin/bash
# Driver for creation of Xenbase gaf

VENV_HOME="$HOME/myenv"    # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# Set to false if file download not required
GET_FILES=false

# -------------------- Download Files -------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps, download xenbase files
if [ "$GET_FILES" = true ]; then
    python3 get_files.py --xen_download --log    # Fresh download of files needed to build xenbase gaf
else
    echo "Skipping file download"
fi

# --------------------- GAF Creation ---------------------

# Parse through xenopus gaf (from EBI), find xenbase matches, & create xenbase gaf
python3 goa_parsing.py --create_gaf --log

# GAF v2.2 post processing
OUTPUT_DIR=$HOME/xenbase-gaf-pipeline/output-files
xb_gaf=$OUTPUT_DIR/xenbase.EBI.only.gaf
./gaf2.2_mods.sh $OUTPUT_DIR $xb_gaf