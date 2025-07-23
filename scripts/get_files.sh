#!/bin/bash
DATE=`date +%Y-%m-%d`

DATA_DIR="$HOME/xenbase-gaf-pipeline/input-files"
GAF_DIR="$DATA_DIR/goa-gafs"
NCBI_MAP_DIR="$DATA_DIR/ncbi-maps"
XB_DIR="$DATA_DIR/xenbase-files"

mkdir -p $GAF_DIR
mkdir -p $GAF_DIR/original-goa
mkdir -p $NCBI_MAP_DIR

# Download goa data for all species
curl -o $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
echo "Full GOA downloaded"

# Extract into species specific files
echo "Starting extraction of species-specific annotations from $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz"
zgrep -a 'taxon:8355\|taxon:8364' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Xenopus.GOA.Extracted.$DATE.gaf
echo "Xenopus GOA Extracted"
zgrep -a 'taxon:9606' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Human.GOA.Extracted.$DATE.gaf
echo "Human GOA Extracted"
zgrep -a 'taxon:10090' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Mouse.GOA.Extracted.$DATE.gaf
echo "Mouse GOA Extracted"
zgrep -a 'taxon:10116' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Rat.GOA.Extracted.$DATE.gaf
echo "Rat GOA Extracted"
zgrep -a 'taxon:9031' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Chicken.GOA.Extracted.$DATE.gaf
echo "Chicken GOA Extracted"
zgrep -a 'taxon:7955' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Zebrafish.GOA.Extracted.$DATE.gaf
echo "Zebrafish GOA Extracted"
zgrep -a 'taxon:7227' $GAF_DIR/original-goa/goa_uniprot_all.gaf.gz > $GAF_DIR/Drosophila.GOA.Extracted.$DATE.gaf
echo "Drosophila GOA Extracted"

# Download Xenbase gpi, genepage to gene ID map, and gene ID to ortholog NCBI id map
curl https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz | gunzip > $XB_DIR/Xenbase.gpi
curl https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping_chd.txt > $XB_DIR/Xenbase_Genepage_To_GeneId.txt
# NOTE: ^ line 11651 may be uncorrectly tabbed ("vma22  .L" instead of vma22.L); fix before using
curl https://xenbase-bio1.ucalgary.ca/cgi-bin/reports/genepage_entrez_orthologs.cgi | sed '/./,$!d' > "$NCBI_MAP_DIR/Xenopus_NCBI_Orthologs.tsv"
{ echo -e "umbrella_id\txtrop\thuman\tmouse\trat\tzebrafish\tchicken\tdrosphila\tworm"; tail -n +2 "$NCBI_MAP_DIR/Xenopus_NCBI_Orthologs.tsv"; } > "$NCBI_MAP_DIR/tmp" && mv "$NCBI_MAP_DIR/tmp" "$NCBI_MAP_DIR/Xenopus_NCBI_Orthologs.tsv"
# NOTE: ^ Changed header in ortholog file to match species names used here: "umbrella_id	xtrop	human	mouse	rat	zebrafish	chicken	drosphila	worm"
echo "Xenbase files downloaded"

# Download uniprot to NCBI id mapping files for each ortholog species
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Human_NCBI_Mapping.tsv
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Mouse_NCBI_Mapping.tsv
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/RAT_10116_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Rat_NCBI_Mapping.tsv
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/CHICK_9031_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Chicken_NCBI_Mapping.tsv
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/DANRE_7955_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Zebrafish_NCBI_Mapping.tsv
curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/DROME_7227_idmapping.dat.gz | gunzip | grep 'GeneID' > $NCBI_MAP_DIR/Drosophila_NCBI_Mapping.tsv
echo "NCBI maps downloaded"

# - - - - - - - - - - - - - - - - - - - - Not Used - - - - - - - - - - - - - - - - - - - -

# Species specific gaf sources (excluding xenopus) (faster downloads, but seem to be missing annotations found in goa_uniprot_all.gaf file):
#curl https://current.geneontology.org/annotations/goa_human.gaf.gz | gunzip > $GAF_DIR/Human.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/mgi.gaf.gz | gunzip > $GAF_DIR/Mouse.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/rgd.gaf.gz | gunzip > $GAF_DIR/Rat.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/goa_chicken.gaf.gz | gunzip > $GAF_DIR/Chicken.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/zfin.gaf.gz | gunzip > $GAF_DIR/Zebrafish.GOA.Currated.$DATE.gaf
#curl https://current.geneontology.org/annotations/fb.gaf.gz | gunzip > $GAF_DIR/Drosophila.GOA.Currated.$DATE.gaf

# Xenbase gaf:
#curl https://download.xenbase.org/xenbase/GenePageReports/xenbase.gaf.gz | gunzip > $XB_DIR/Xenbase.gaf

# See all uniport to NCBI mapping files available
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
