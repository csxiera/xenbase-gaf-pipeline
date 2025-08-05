# Xenbase GAF Pipeline
This repository contains scripts for building Gene Ontology Annotation files for Xenopus. Annotations for both Xenopus laevis and Xenopus tropicalis are sourced from the GO Central repository and linked to Xenbase using gene entity identifiers (XB-GENE-IDs). Ortholog annotations are also identified and included in GAF output as 'ISO' annoations. Ortholog species include: human, mouse, rat, chicken, zebrafish, and drosophila

## Scripts
1. `Xenbase_gaf_driver.sh`    - DRIVER: script to run GAF pipeline
2. `get_files.py`             - MAIN: script to download various files used in pipeline
3. `goa_parsing.py`           - MAIN: script used to parse through GOA GAF files, create Xenbase GAF, and find ortholog annotations
4. `compare_files.py`         - HELPER: script to compare 2 files, writes differing lines to seperate files for analysis
5. `count_annotations.py`     - HELPER: script to count the number of annotations in a file, with or without filtering applied

More detailed information about each program is found in the script itself

## Running Pipeline
Xenbase_gaf_driver.sh should be set up to use the current date, with GET_FILES = true when used in a cron job  
- To rerun GAF pipeline without redownloading files, change current date to previous download date and set GET_FILES = false  

Input folder structure is set up in `get_files.py`. Output folder structure is set up in `goa_parsing.py`. Program expects `xenbase-gaf-pipeline` project folder to be in HOME directory  

Helper scripts are designed to be run as standalone programs. The `__main__` section contains example calls demonstrating how to call the helper functions for different usages.

## Input/Output

### Required inputs:
1. GOAs:
    - {species} goa extracted GAF (x7; one for each species)
2. NCBI to UniProt ID Maps:
    - {species} NCBI mapping (x7, one for each species)
3. Xenbase files:
    - Xenbase GPI
    - Xenbase genepage to gene ID map

### Expected outputs:
1. GAFs:
    - Xenbase
    - Xenbase (x.trop only)
    - Xenbase (x.laev only)
    - Xenbase with all orthologs
    - Xenbase with {species} (human/mouse/rat/chicken/zebrafish/drosophila)
    - All orthologs
2. TSVs:
    - Unmatched annotations
    - Annotation provenance

Outputs from helper scripts can be found in `xenbase/output-files/analysis`

---
### Sources

Adapted from original Xenbase GAF pipeline:
https://gitlab.com/Xenbase/GOA_pipeline.git

GOA GAF annotations for each species extracted from:
https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

NCBI to UniProt Mapping files from:
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/