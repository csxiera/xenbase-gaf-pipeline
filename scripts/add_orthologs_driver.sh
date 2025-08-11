#!/bin/bash
# Driver for adding ortholog annotatiosn to post-NOCTUA xenbase gaf

# Parse through ortholog gafs (from EBI), find matches to xenbase genes, add as ISO annotations to xenbase GAF
python3 goa_parsing.py --date --add_orthos $DATE #--log