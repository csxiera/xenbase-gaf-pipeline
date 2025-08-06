#!/usr/bin/perl
# From original GOA pipeline (see 'scripts/original-goa-scripts')
# Modified for output comparison with new GAF processing pipeline

# use strict;
# use warnings;
use File::Basename;

my $GOA=$ARGV[0];   # Define GOA input file as first command line argument
my $DATE=$ARGV[1];
my $goa_basename = basename($GOA);
my $base_dir = "$ENV{HOME}/xenbase-gaf-pipeline";
my $gpiungz = "$base_dir/input-files/xenbase-files/Xenbase.gpi";
my $out_dir = "$base_dir/output-files/original-pipeline-output";

print "\tMapping UniProts to Xenbase from GPI...\n";
open(GPI, "<:encoding(UTF-8)", $gpiungz) or die "Could not open $gpiungz: $!";
open (GOA, "<:encoding(UTF-8)", $GOA);#Read in current GOA file
while (<GPI>){
	next if $_=~m/^!/;#Ignore commented out lines in GPI
	chomp;
	my($xenbase,$DB_symbol,$DB_Object_Name,$DB_Object_Synonym,$dbxref)=(split("\t"))[1,2,3,4,8];#Pull Xenbase ID, object name, synonym and DBxref columns from GPI. The current setting of column 8 for the DBxrefs is for the v1.2 gpi format, if we transition to v2 the final column pulled will need changed to 10.
	my(@uniprots)=($dbxref=~m/Uni[Pp]rotKB:(.+?)(?:\||$)/gi);#find all 'UniProtKB:######' instances and make into array.
	# my (@symbols)=($DB_symbol);# make array for symbols
	# my (@Gene_names)=($DB_symbol);# make array for gene names
	# my (@Gene_synonyms)=($DB_Object_Synonym);# make array for synonyms
	foreach $uniprot(@uniprots){
		$uniprot2Xenbase{$uniprot}=$xenbase;
	}
	foreach $uniprot(@uniprots){
		$uniprot2symbols{$uniprot}=$DB_symbol;
	}
		foreach $uniprot(@uniprots){
		$uniprot2gene_name{$uniprot}=$DB_Object_Name;
	}
	foreach $uniprot(@uniprots){
		$uniprot2gene_synonyms{$uniprot}=$DB_Object_Synonym;
	}
}
print "\tFinished mapping.\n";

print "\n\tMatching GAF annotations...\n";
open (my $matched,">:encoding(UTF-8)","$out_dir/matched_annotations.Xenopus.GOA.$DATE.gaf") or die $!;
print $matched "!gaf-version: 2.2
!Project_name: Xenbase 
!URL: http://www.xenbase.org/
!Contact Email: xenbase\@cchmc.org
!Funding: NICHD of US National Institutes of Health
!date: $DATE 
!
!	Column	Content	Required?	Cardinality	Example
!	1	DB	required	1	Xenbase
!	2	DB Object ID	required	1	XB-GENE-920512
!	3	DB Object Symbol	required	1	ventx1.1.L
!	4	Qualifier	required	1 or 2	NOT|acts_upstream_of_or_within
!	5	GO ID	required	1	GO:0030509
!	6	DB:Reference (|DB:Reference)	required	1 or greater	PMID:11944926
!	7	Evidence Code	required	1	IMP
!	8	With (or) From	optional	0 or greater	GO:0043565
!	9	Aspect	required	1	P
!	10	DB Object Name	optional	0 or 1	VENT homeobox 1, gene 1 L homeolog
!	11	DB Object Synonym (|Synonym)	optional	0 or greater	PV.1|Xvent-1b|posterior-ventral 1|ventx1.1-a|ventx1.1-b
!	12	DB Object Type	required	1	Protein
!	13	Taxon(|taxon)	required	1 or 2	taxon:8355
!	14	Date	required	1	20090118
!	15	Assigned By	required	1	Xenbase
!	16	Annotation Extension	optional	0 or greater	part_of(CL:0000576)
!	17	Gene Product Form ID	optional	0 or 1	UniProtKB:P12345-2
!\n";
close $matched;

while (<GOA>) {
	next if $_=~m/^!/; # Ignore commented out lines in GOA
	chomp;
	# my($DB,$UNIPROTID,$symbol,$qualifier,$GOID,$DBref,$evidence,$withfrom,$other_annotation)=(split("\t",$_,9))[0,1,2,3,4,5,6,7,8];#Older section just for substituting XB-GENE-IDs for Uniprot IDs.
	my($DB,$UNIPROTID,$symbol,$qualifier,$GOID,$DBref,$evidence,$withfrom,$Aspect,$object_name,$object_synonyms,$object_type,$taxon,$date,$assigned_by,$annotation_extension,$gene_product_id)=(split("\t",$_,17))[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];#modified for changing more elements in GAF
	my $XBGENE=$uniprot2Xenbase{$UNIPROTID};
	my $XB_symbol=$uniprot2symbols{$UNIPROTID};
	my $XB_gene_name=$uniprot2gene_name{$UNIPROTID};
	my $XB_gene_synonyms=$uniprot2gene_synonyms{$UNIPROTID};
	if ($XBGENE=~m/XB-GENE/) { # This lines only prints out results if there is a matched XB-GENE-ID in column two.
		open (my $matched,">>:encoding(UTF-8)","$out_dir/matched_annotations.Xenopus.GOA.$DATE.gaf");
		print $matched "Xenbase\t$XBGENE\t$XB_symbol\t$qualifier\t$GOID\t$DBref\t$evidence\t$withfrom\t$Aspect\t$XB_gene_name\t$XB_gene_synonyms\t$object_type\t$taxon\t$date\t$assigned_by\t$annotation_extension\t$gene_product_id\n";
		open (my $provenance,">>:encoding(UTF-8)","$out_dir/annotation_provenance.tmp");
		print $provenance "Xenbase\t$XBGENE\t$XB_symbol\t$qualifier\t$GOID\t$DBref\t$evidence\t$withfrom\t$Aspect\t$XB_gene_name\t$XB_gene_synonyms\t$object_type\t$taxon\t$date\t$assigned_by\t$annotation_extension\t$DB:$UNIPROTID\n";# This elsif section checks if the withfrom field is filled
		}
	else{
        open (my $unmatched,">>:encoding(UTF-8)", "$out_dir/unmatched_annotations.Xenopus.GOA.$DATE.tmp");
        print $unmatched "$DB\t$UNIPROTID\t$symbol\t$qualifier\t$GOID\t$DBref\t$evidence\t$withfrom\t$Aspect\t$object_name\t$object_synonyms\t$object_type\t$taxon\t$date\t$assigned_by\t$annotation_extension\t$gene_product_id\n";
	}
}
print "\tFinished matching annotations.\n";