#!/bin/bash
DATE=`date +%Y-%m-%d`

DATA_DIR=~/xenbase-gaf-pipeline/input-files
GAF_DIR="$DATA_DIR/gaf"
NCBI_MAP_DIR="$DATA_DIR/ncbi-map"

mkdir -p $GAF_DIR
mkdir -p $NCBI_MAP_DIR

# Download goa data for all species, extract into species specific files
#curl https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz |gunzip -c|grep -a 'taxon:8355\|taxon:8364' > $GAF_DIR/Xenopus.GOA.Extracted.$DATE.gaf
#curl https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz |gunzip -c|grep -a 'taxon:9606' > $GAF_DIR/Human.GOA.Extracted.$DATE.gaf

zgrep -a 'taxon:8355\|taxon:8364' https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz > $GAF_DIR/Xenopus.GOA.Extracted2.$DATE.gaf

# Species specific gaf sources (excluding xenopus) (faster downloads, but seem to be missing annotations found in goa_uniprot_all.gaf file):
#curl https://current.geneontology.org/annotations/goa_human.gaf.gz | gunzip > $GAF_DIR/Human.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/mgi.gaf.gz | gunzip > $GAF_DIR/Mouse.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/rgd.gaf.gz | gunzip > $GAF_DIR/Rat.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/goa_chicken.gaf.gz | gunzip > $GAF_DIR/Chicken.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/zfin.gaf.gz | gunzip > $GAF_DIR/Zebrafish.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/fb.gaf.gz | gunzip > $GAF_DIR/Drosophila.GOA.Currated.$DATE.gaf

# Download Xenbase gpi, genepage to gene ID map, and gene ID to ortholog NCBI id map
#curl https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz | gunzip > $DATA_DIR/Xenbase.gpi
#curl https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping_chd.txt > $DATA_DIR/Xenbase_Genepage_To_GeneId.txt
#curl https://xenbase-bio1.ucalgary.ca/cgi-bin/reports/genepage_entrez_orthologs.cgi > $NCBI_MAP_DIR/Xenopus_NCBI_Orthologs.tsv

# Download uniprot to NCBI id mapping files for each ortholog species
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Human_NCBI_Mapping.tsv
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Mouse_NCBI_Mapping.tsv




# See all uniport to NCBI mapping files available
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
