#!/bin/bash
# Driver for creation of Xenbase gaf

VENV_HOME="$HOME/myenv"    # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# -------------------- Set Variables --------------------
# Set date to use: 
#DATE=$(date +"%Y-%m-%d")    # Current date for new downloads
DATE='2025-07-18'          # Past date to use previously downloaded files

# Set to true if download or extraction is required
GET_FILES=false

# -------------------- Download Files -------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps, download xenbase files
if [ "$GET_FILES" = true ]; then
    TODAY=$(date +"%Y-%m-%d")

    if [ "$DATE" = "$TODAY" ]; then
        python3 get_files.py --xen_goa --ortho_goa --maps --xenbase --log    # Fresh download of all required files
    else
        python3 get_files.py --extract --date "$DATE" --log                  # Re-extraction from previous goa download (assumes mapping/xenbase files already downloaded)
    fi
else
    echo "Skipping file download and extraction"
fi

# --------------------- GAF Creation ---------------------

# Parse through xenopus gaf (from EBI), find xenbase matches, & create xenbase gaf
python3 goa_parsing.py --date --create_gaf $DATE #--log

# GAF v2.2 post processing
OUTPUT_DIR=$HOME/xenbase-gaf-pipeline/output-files
xb_gaf=$OUTPUT_DIR/xenbase.EBI.only.gaf
./gaf2.2_mods.sh $OUTPUT_DIR $xb_gaf