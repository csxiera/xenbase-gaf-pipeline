# Xenbase GAF Pipeline
This repository contains scripts for building Gene Ontology Annotation files for Xenopus. Annotations for both Xenopus laevis and Xenopus tropicalis are sourced from the GO Central repository and linked to Xenbase using gene entity identifiers (XB-GENE-IDs). Ortholog annotations are also identified and included in GAF output as 'ISO' annoations. Ortholog species include: human, mouse, rat, chicken, zebrafish, and drosophila

## Scripts
1. `.gitlab-ci.yml`           - DRIVER: script to run full pipeline
2. `xenbase_gaf_driver.sh`    - DRIVER: script to build Xenbase GAF from GPI & GOA files
3. `add_orthologs_driver.sh`  - DRIVER: script to add ortholog annotations to Xenbase GAF
4. `get_files.py`             - MAIN: script to download various files used in pipeline
5. `goa_parsing.py`           - MAIN: script used to parse through GOA GAF files, create Xenbase GAF, and find ortholog annotations
6. `gaf2.2_mods.sh`           - MAIN: script used to modify Xenbase GAF output to be consistent with GAF version 2.2
7. `compare_files.py`         - HELPER: script to compare 2 files, writes differing lines to seperate files for analysis
8. `count_annotations.py`     - HELPER: script to count the number of annotations in a file, with or without filtering applied

NOTE: The subfolder `original-pipeline` contains adapted version of the original GOA Pipeline bash/perl scripts. These were used to verify that the `xenbase.EBI.only.2.2.gaf` output from both pipelines are identical.

More detailed information about each program can be found in the scripts themselves. Project flow can be viewed in `Pipeline Flowchat.png`

## Running Pipeline
Pipeline is designed to run via CI/CD pipeline, which will:  
1. Build Xenbase GAF on scheduled or manual trigger
2. Add orthologs to the post-NOCTUA Xenbase GAF when this file is added/changed

Both `xenbase_gaf_driver.sh` and `add_orthologs_driver.sh` should be set up with `GET_FILES = true`
- For testing purposes: to rerun GAF pipeline without redownloading files, set `GET_FILES = false ` 

Input folder structure is set up in `get_files.py`. Output folder structure is set up in `goa_parsing.py`. Program expects `xenbase-gaf-pipeline` repo folder to be in `HOME` directory  

Helper scripts are designed to be run as standalone programs for testing purposes. The `__main__` section contains example calls demonstrating how to call the helper functions for different usages.

## Input/Output

### Required inputs:
1. GOAs:
    - x.trop goa GAF
    - x.laev goa GAF
    - {species} goa extracted GAF (x6; one for each ortholog species)
2. NCBI to UniProt ID Maps:
    - {species} NCBI mapping (x7, one for each species)
3. Xenbase files:
    - Xenbase GPI
    - Xenbase genepage to gene ID map

### Expected outputs:
1. GAFs:
    - Xenbase
    - Xenbase x.trop only
    - Xenbase x.laev only
    - Xenbase plus all orthologs
    - Xenbase plus {species} (human/mouse/rat/chicken/zebrafish/drosophila)
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