# ------------------- GAF 2.2 Formatting -------------------

OUTPUT_DIR=$1
input_gaf=$2
output_gaf="${input_gaf%.gaf}.2.2.gaf"

echo -e "Formatting for GAF v2.2...\n"

# Filters out GOC and GO_Central annotations
cat $input_gaf | grep -av '^!' |grep -av 'GO_Central'|grep -av 'GOC'| sort | uniq > $OUTPUT_DIR/non-redundant.gaf.tmp

# Steps 1-4 make the qualifier columns consistent with default GAF 2.2 relationships and changes the targets of self binding to pass GO quality checks:
tmp=$OUTPUT_DIR/non-redundant.gaf.tmp

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

# Write header + deduplicate
grep -a '^!' "$input_gaf" > "$output_gaf"
sort "$tmp" | uniq >> "$output_gaf"

# Clean up
rm "$tmp"
rm "$input_gaf"