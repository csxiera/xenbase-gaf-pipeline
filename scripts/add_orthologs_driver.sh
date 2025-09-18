#!/bin/bash
# Driver for adding ortholog annotations to post-NOCTUA xenbase gaf

VENV_HOME="$HOME/myenv"     # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# Set ortho species to process
ortho_species='["Human", "Mouse", "Rat", "Chicken", "Zebrafish", "Drosophila"]'     # All ortholog species currently supported

# Set to false if file download/extraction not required
GET_FILES=true

# -------------------- Download Files -------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps
if [ "$GET_FILES" = true ]; then
        python3 get_files.py --ortho_download --log         # Fresh download of files needed to add orthologs
else
    echo "Skipping file download and extraction"
fi

# ------------------- GAF Modification -------------------

# Parse through ortholog gafs (from EBI), find matches to xenbase genes, add as ISO annotations to xenbase GAF
python3 goa_parsing.py --add_orthos "$ortho_species" --log