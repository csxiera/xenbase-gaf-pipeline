#!/bin/bash
# From original GOA pipeline
# Modified for output comparison with new GAF processing pipeline

MAIN_DIR=$HOME/xenbase-gaf-pipeline
SCRIPT_DIR=$MAIN_DIR/scripts/original-pipeline
INPUT_DIR=$MAIN_DIR/input-files
OUTPUT_DIR=$MAIN_DIR/output-files/original-pipeline-output
DATE='2025-07-18'   # Download/extraction date of GOA files to use

XB_GAF=$INPUT_DIR/goa-gafs/Xenopus.GOA.Extracted.$DATE.gaf

echo -e "\nRunning parsing script...\n"
perl $SCRIPT_DIR/GOA_parsing_adapted.pl $XB_GAF $DATE
echo -e "\nParsing finished.\n"

echo -e "Removing redundant lines...\n"
matched=$OUTPUT_DIR/matched_annotations.Xenopus.GOA.$DATE.gaf
cat $matched | grep -av '^!' |grep -av 'GO_Central'|grep -av 'GOC'| sort | uniq > $OUTPUT_DIR/non-redundant.combined.Xenopus.$DATE.gaf.tmp # This line filters out GOC and GO_Central annotations.

# Steps 1-4 make the qualifier columns consistent with default gaf2.2 relationships and changes the targets of self binding to pass GO quality checks:
echo -e "Formatting for GAF v2.2...\n"
non_redundant=$OUTPUT_DIR/non-redundant.combined.Xenopus.$DATE.gaf.tmp
tmp=$OUTPUT_DIR/non-redundant.Xenopus.tmp
cp $non_redundant $tmp

# Step 1a: If qualifier (column 4) is 'involved_in' and aspect (column 9) is 'P', change qualifer to 'acts_upstream_of_or_within'
awk -F'\t' 'BEGIN{OFS="\t"}{if($4=="involved_in" && $9=="P"){$4="acts_upstream_of_or_within"}{print}}' $tmp > ${tmp}.new && mv ${tmp}.new $tmp
# Step 1b: Same as above but for 'NOT|' qualifier
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
sort $tmp | uniq > $OUTPUT_DIR/Xenbase.EBI.only.$DATE.gaf
rm $non_redundant $tmp

# Copy GAF header from matched file to EBI only file
ebi_final=$OUTPUT_DIR/Xenbase.EBI.only.2.2.gaf
cat $matched|grep -a '^!' | cat - $OUTPUT_DIR/Xenbase.EBI.only.$DATE.gaf > $ebi_final
rm $OUTPUT_DIR/Xenbase.EBI.only.$DATE.gaf

# Get summary stats for EBI GAF final file
ebi_final_summary=$OUTPUT_DIR/EBI_only_2.2_Summary.txt
echo -e 'Contributions by source:'\\n"$( cat $ebi_final|grep -av '\!'| cut -f15| sort| uniq -c|sed 's/^[ ]*//g'|awk '{print $2"\t"$1}' )" > $ebi_final_summary # contributors
echo -e 'Contributions by Type:'\\n"$( cat $ebi_final|grep -av '\!'| cut -f7| sort| uniq -c|sed 's/^[ ]*//g'|awk '{print $2"\t"$1}' )" >> $ebi_final_summary # evidence codes
echo -e 'Unique Genes: ' $( cat $ebi_final|grep -av '\!'| cut -f2| sort| uniq |wc -l) >> $ebi_final_summary # Unique genes
echo -e 'Unique GO terms: '$( cat $ebi_final|grep -av '\!'| cut -f5| sort| uniq |wc -l ) >> $ebi_final_summary # Unique GO terms
echo -e 'Unique Gene-GO term associations: '$( cat $ebi_final|grep -av '\!'| cut -f2,5| sort| uniq |wc -l ) >> $ebi_final_summary #Unique gene GO associations
filesize=$( ls -lSr $ebi_final | awk '{print $5}' | head -n1 ) # The appropriate column may differ depending on system, this works on gitlab runner with perl docker.
echo $ebi_final file is $filesize bytes
error=$( expr $filesize \>=  2000 ) # Checks if any of the GAFs have a size of less than 2Kb and fails the job if they do
if [ "$error" == "0" ]
then 
    echo 'GAF file is a stub' 
	exit 1
else 
	echo 'GAF file seems full length'
fi

echo -e "\nFinished!"