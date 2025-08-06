#!/bin/bash
# Driver for creation of Xenbase gaf and Xenbase with ortholog gafs

VENV_HOME="$HOME/myenv"    # Modify path to point to virtual environment directory
source $VENV_HOME/bin/activate

# -------------------- Set Variables --------------------
# Set date to use: 
#DATE=$(date +"%Y-%m-%d")    # Current date for new downloads
DATE='2025-07-18'          # Past date to use previously downloaded files

# Set to true if download or extraction is required
GET_FILES=false
# --------------------- GAF Creation ---------------------

# Download goa gaf, extract to species specific gafs, download ncbi -> uniprot maps, download xenbase files
if [ "$GET_FILES" = true ]; then
    TODAY=$(date +"%Y-%m-%d")

    if [ "$DATE" = "$TODAY" ]; then
        python3 get_files.py --goa --extract --maps --xenbase --log    # Fresh download of all required files
    else
        python3 get_files.py --extract --date "$DATE" --log            # Re-extraction from previous goa download (assumes mapping/xenbase files already downloaded)
    fi
else
    echo "Skipping file download and extraction"
fi

# Parse through gaf files, find xenbase matches, find ortholog annotations, combine xenbase & ortholog annotations
python3 goa_parsing.py --date $DATE #--log

# ------------------- GAF 2.2 Formatting -------------------

OUTPUT_DIR=$HOME/xenbase-gaf-pipeline/output-files
xb_gaf=$OUTPUT_DIR/Xenbase.gaf

echo -e "Formatting for GAF v2.2...\n"

# Filters out GOC and GO_Central annotations
cat $xb_gaf | grep -av '^!' |grep -av 'GO_Central'|grep -av 'GOC'| sort | uniq > $OUTPUT_DIR/non-redundant.Xenbase.gaf.tmp

# Steps 1-4 make the qualifier columns consistent with default GAF 2.2 relationships and changes the targets of self binding to pass GO quality checks:
tmp=$OUTPUT_DIR/non-redundant.Xenbase.gaf.tmp

# Step 1a: If qualifier (column 4) is 'involved_in' and aspect (column 9) is 'P', change qualifer to 'acts_upstream_of_or_within'
awk -F'\t' 'BEGIN{OFS="\t"}{if($4=="involved_in" && $9=="P"){$4="acts_upstream_of_or_within"}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp
# Step 1b: Same as above but for 'NOT|' qualifiers
awk -F'\t' 'BEGIN{OFS="\t"}{if($4=="NOT|involved_in" && $9=="P"){$4="NOT|acts_upstream_of_or_within"}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp

# Step 2a: Adds 'taxon:' prefix to taxon ID 8355 if missing
awk -F'\t' 'BEGIN{OFS="\t"}{if($13=="8355"){$13="taxon:8355"}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp
# Step 2b: Same as above but for taxon ID 8364
awk -F'\t' 'BEGIN{OFS="\t"}{if($13=="8364"){$13="taxon:8364"}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp

# Step 3a: If the GO term (column 5) is 'GO:0042803' and the evidence code (column 7) IS NOT in [ISS/ISA/ISO/IDA], set with/from (column 8) to 'Xenbase:{gene-ID}'
awk -F'\t' 'BEGIN{OFS="\t"}{if($5=="GO:0042803" && !($7=="ISS"||$7=="ISA"||$7=="ISO"||$7=="IDA")){$8="Xenbase:"$2}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp
# Step 3b: Same as above but for 'GO:0051260'
awk -F'\t' 'BEGIN{OFS="\t"}{if($5=="GO:0051260" && !($7=="ISS"||$7=="ISA"||$7=="ISO")){$8="Xenbase:"$2}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp

# Step 4a: If the GO term is 'GO:0042803' and the evidence code IS in [ISS/ISA/ISO], then append 'Xenbase:{gene-ID}|' to existing with/from
awk -F'\t' 'BEGIN{OFS="\t"}{if($5=="GO:0042803" && ($7=="ISS"||$7=="ISA"||$7=="ISO")){$8="Xenbase:"$2"|" $8}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp
# Step 4b: Same as above but for 'GO:0051260'
awk -F'\t' 'BEGIN{OFS="\t"}{if($5=="GO:0051260" && ($7=="ISS"||$7=="ISA"||$7=="ISO")){$8="Xenbase:"$2"|" $8}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp

# Sort, deduplicate, and remove temp files
sort $tmp | uniq > $OUTPUT_DIR/Xenbase.EBI.only.gaf
rm $tmp

# Copy GAF header from matched file to EBI only file
ebi_final=$OUTPUT_DIR/Xenbase.EBI.only.2.2.gaf
cat $xb_gaf|grep -a '^!' | cat - $OUTPUT_DIR/Xenbase.EBI.only.gaf > $ebi_final
rm $OUTPUT_DIR/Xenbase.EBI.only.gaf