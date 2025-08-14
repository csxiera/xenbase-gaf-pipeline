#!/bin/bash
# Driver for creation of Xenbase gaf

VENV_HOME="$HOME/myenv"    # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# -------------------- Set Variables --------------------
# Set date to use: 
DATE=$(date +"%Y-%m-%d")    # Current date for new downloads
#DATE='2025-07-18'          # Past date to use previously downloaded files

# Set to false if file download not required
GET_FILES=true

# -------------------- Download Files -------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps, download xenbase files
if [ "$GET_FILES" = true ]; then
    python3 get_files.py --xen_download --date $DATE --log    # Fresh download of files needed to build xenbase gaf
else
    echo "Skipping file download"
fi

# --------------------- GAF Creation ---------------------

# Parse through xenopus gaf (from EBI), find xenbase matches, & create xenbase gaf
python3 goa_parsing.py --create_gaf --date $DATE --log

# GAF v2.2 post processing
OUTPUT_DIR=$HOME/xenbase-gaf-pipeline/output-files
xb_gaf=$OUTPUT_DIR/xenbase.EBI.only.gaf
./gaf2.2_mods.sh $OUTPUT_DIR $xb_gaf