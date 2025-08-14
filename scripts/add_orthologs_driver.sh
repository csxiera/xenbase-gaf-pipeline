#!/bin/bash
# Driver for adding ortholog annotatiosn to post-NOCTUA xenbase gaf

VENV_HOME="$HOME/myenv"    # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# -------------------- Set Variables --------------------
# Set date to use: 
DATE=$(date +"%Y-%m-%d")    # Current date for new downloads
#DATE='2025-07-18'          # Past date to use previously downloaded files

# Set ortho species to process
ortho_species='["Human", "Mouse", "Rat", "Chicken", "Zebrafish", "Drosophila"]'     # All ortholog species currently supported

# Set to false if file download or extraction not required
GET_FILES=true

# -------------------- Download Files -------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps
if [ "$GET_FILES" = true ]; then
    TODAY=$(date +"%Y-%m-%d")

    if [ "$DATE" = "$TODAY" ]; then
        python3 get_files.py --ortho_download --date $DATE --log           # Fresh download of files needed to add orthologs
    else
        python3 get_files.py --extract --date "$DATE" --log                # Re-extraction from previous goa download (assumes mapping/xenbase files already downloaded)
    fi
else
    echo "Skipping file download and extraction"
fi

# ------------------- GAF Modification -------------------

# Parse through ortholog gafs (from EBI), find matches to xenbase genes, add as ISO annotations to xenbase GAF
python3 goa_parsing.py --add_orthos "$ortho_species" --date $DATE --log