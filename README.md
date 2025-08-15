# Xenbase GAF Pipeline
This repository contains scripts for building Gene Ontology Annotation files for Xenopus. Annotations for both Xenopus laevis and Xenopus tropicalis are sourced from the GO Central repository and linked to Xenbase using gene entity identifiers (XB-GENE-IDs). Ortholog annotations are also identified and included in GAF output as 'ISO' annoations. Ortholog species include: human, mouse, rat, chicken, zebrafish, and drosophila

## Scripts
1. `run_pipeline.yaml`        - DRIVER: script to run full pipeline
1. `xenbase_gaf_driver.sh`    - DRIVER: script to build Xenbase GAF from GPI & GOA files
2. `add_orthologs_driver.sh`  - DRIVER: script to add ortholog annotations to Xenbase GAF
3. `get_files.py`             - MAIN: script to download various files used in pipeline
4. `goa_parsing.py`           - MAIN: script used to parse through GOA GAF files, create Xenbase GAF, and find ortholog annotations
5. `gaf2.2_mods.sh`           - MAIN: script used to modify Xenbase GAF output to be consistent with GAF version 2.2
6. `compare_files.py`         - HELPER: script to compare 2 files, writes differing lines to seperate files for analysis
7. `count_annotations.py`     - HELPER: script to count the number of annotations in a file, with or without filtering applied

NOTE: The subfolder `original-pipeline` contains adapted version of the original GOA Pipeline bash/perl scripts. These were used to verify whether the `xenbase.gaf` outputs from both pipelines are identical when provided with the same input files.  

More detailed information about each program is found in the script itself. Project flow can be viewed in `Pipeline Flowchat.png`

## Running Pipeline
Pipeline is designed to run via `run_pipeline.yaml`, which will:  
1. Build Xenbase GAF on scheduled or manual trigger
2. Add orthologs to the post-NOCTUA Xenbase gaf once pushed to pipeline

The drivers called in the yaml script download different files:
- `xenbase_gaf_driver.sh` downloads both X.trop and X.laev GOA files and the Xenbase GPI file
- `add_orthologs_driver.sh` downloads the full GOA file (& extracts into species-specific GOA files), all ncbi mapping files, and the Xenbase genepage to gene id mapping file  
NOTE: The Xenbase GPI file is required to add ortholog annotations, therefore it is assumed that `xenbase_gaf_driver.sh` has already run and downloaded that file

Both `xenbase_gaf_driver.sh` and `add_orthologs_driver.sh` should be set up with `GET_FILES = true`
- For testing purposes: to rerun GAF pipeline without redownloading files, set `GET_FILES = false ` 

Input folder structure is set up in `get_files.py`. Output folder structure is set up in `goa_parsing.py`. Program expects `xenbase-gaf-pipeline` repo folder to be in `HOME` directory  

Helper scripts are designed to be run as standalone programs for testing purposes. The `__main__` section contains example calls demonstrating how to call the helper functions for different usages.

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