#!/bin/bash
DATE=`date +%Y-%m-%d`

mkdir -p ~/tmp
cd ~/tmp

# Download goa data for all species
#curl 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
gunzip -c goa_uniprot_all.gaf.gz > goa_uniprot_all.$DATE.gaf

# Seperate gaf for each species
grep -a 'taxon:8355\|taxon:8364' goa_uniprot_all.$DATE.gaf > Xenopus.GOA.$DATE.gaf
grep -a 'taxon:9606' goa_uniprot_all.$DATE.gaf > Human.GOA.$DATE.gaf

# Download Xenbase gpi, gaf, genepage to gene ID map, and gene ID to ortholog NCBI id map
#curl https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz | gunzip > xenbase.gpi
#curl https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping_chd.txt
#curl https://download.xenbase.org/xenbase/GenePageReports/xenbase.gaf.gz | gunzip > xenbase.gaf
#curl https://xenbase-bio1.ucalgary.ca/cgi-bin/reports/genepage_entrez_orthologs.cgi > xenopus_orthologs.tsv

# Download uniprot to NCBI id mapping files for each ortholog species
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz | gunzip | grep 'GeneID' > HUMAN_9606_idmapping.NCBI.tsv
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz | gunzip | grep 'GeneID' > MOUSE_10090_idmapping.NCBI.tsv

# ---------------------- Not used/needed in goa_parsing.py program ----------------------

# Gaf source for colab notebook scripts:
#curl https://current.geneontology.org/annotations/goa_human.gaf.gz

# See all uniport to NCBI mapping files available
#curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
