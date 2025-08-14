import argparse
import ast
import os
import sys
import gzip
import shutil
import requests
import re
import csv
import textwrap
from datetime import date, datetime

# AUTHOR: C. Lenz
#
# SCRIPT FUNCTION: 
#   1. Create maps of uniprot <-> xenbase, xenbase <-> ortholog NCBIs, and uniprot <-> ortholog NCBIs (when needed)
#   2. Create Xenbase gaf from xenopus goa gaf:
#       a. filter out uniprots not mapped to Xenbase gene IDs
#       b. replace uniprot ID with Xenbase gene IDs
#       c. record unmatched annotations
#   3. For each ortholog species gaf:
#       a. filter out non-evidence based codes 
#       b. filter out annotations with uniprots not mapped to xenbase
#       c. replace ortholog uniprot ID with Xenbase gene ID (for X.trop) move ortholog uniprot ID to 'with/from' field
#   4. Collapse ortholog annotations into unique gene + GO term pairs, with ALL ortholog uniprots listed in 'with/from' field
#   5. Add collapsed annotations to Xenbase gaf (after processing by NOCTUA pipeline)
#
# NOTE: Original GOA_parsing.pl script translated & adapted into parse_gaf, create_xenbase_gaf, and match_to_xenbase functions

# ---------------------------------------- Main Functions ----------------------------------------

# FUNCTION: Main function to create xenbase GAF
def main_gaf(dl_date):
    # Date xenopus files were downloaded/extracted
    print(f"Download date of xenopus GOA file used: {dl_date}")

    # Populate Xenopus maps
    xenopus_gaf = os.path.join(gaf_dir, f'Xenopus.GOA.Curated.{dl_date}.gaf')
    populate_maps(xenopus_gaf)

    # Create Xenbase GAF from GPI & GOA
    create_xenbase_gaf(xenopus_gaf, output_dir)

# FUNCTION: Main function to add ortholog 'ISO' annotations into xenbase GAF
def main_ortho(dl_date, ortho_species):
    # Date ortholog GOA files were downloaded/extracted
    print(f"Download date of ortholog GOA files used: {dl_date}")

    # Split post-noctua xenbase gaf into x.trop & x.laev files                  # CHECK!!: noctua output may need unzipping first
    xenbase_gaf = os.path.join(output_dir, "xenbase.EBI.only.2.2.gaf")          # FIX!!: Change to noctua output filename if different
    populate_maps(xenbase_gaf, split=True)
    split_xenbase_gaf(xenbase_gaf, output_dir, zip=True)

    # Populate xenbase <-> ortholog lookup maps
    populate_maps(xenbase_gaf, ortho_lookup=True)
    os.makedirs(os.path.join(output_dir, "ortho-gafs"), exist_ok=True)

    # Modify ortholog GAFs for each species with xenbase IDs
    for species in ortho_species:
        # Clear ortholog maps (for clean population)
        ortho_uniprot_map.clear()

        print(f"\n--------------- Processing {species} Ortholog Annotations ---------------")

        # Filter ortholog GAF files to remove invalid evidence codes before populating maps
        gaf_suffix = f'GOA.Extracted.{dl_date}.gaf'
        ortho_gaf = os.path.join(gaf_dir, f'{species}.{gaf_suffix}')
        filter_gaf(ortho_gaf, create_copy=True, filter_col=6, filter_values=allowed_codes)        

        # Populate ortholog maps for given species
        filtered_ortho_gaf = os.path.join(gaf_dir, f'{species}.{gaf_suffix}')
        populate_maps(filtered_ortho_gaf, species=species)

        # Filter ortholog gaf to only include annotations with xenopus ncbi gene id matches
        filtered_ortho_gaf = filter_gaf(filtered_ortho_gaf, create_copy=False, filter_col=1, filter_values=ortho_uniprot_map.keys())
        
        # Map ortholog annotations to Xenbase gene IDs and sub Xenbase info into ortholog GAF
        xen_from_ortho = modify_ortho_annotations(filtered_ortho_gaf, output_dir, species)

        # Combine ortholog annotations with Xenbase GAF
        xen_w_ortho = os.path.join(output_dir, f"xenbase.plus.{species.lower()}.gaf")
        combine_annotations(xenbase_gaf, xen_from_ortho, xen_w_ortho)

        # Compress to .gz, remove uncompressed version
        gzip_file(xen_w_ortho, delete_input=True)

    # Remove filtered gafs (.tmp) from input goa gaf folder
    remove_tmp_files(gaf_dir)

    # Create master ortholog file to compile all ortholog annotations in
    master_ortho_gaf = f'{output_dir}/ortho-gafs/Master_Orthologs.gaf'
    open(master_ortho_gaf, 'w').close()

    print(f"\n--------------- Adding Ortholog Annotations ---------------\n")

    for species in ortho_species:
        print(f"Adding {species} ISO annotations to master ortholog file...\n")
        single_ortho_gaf = f'{output_dir}/ortho-gafs/Xenbase_from_{species}.gaf'

        # Add annotations from given species into master ortholog file
        combine_annotations(master_ortho_gaf, single_ortho_gaf, master_ortho_gaf)

        # Collapse duplicate lines by combining uniprots/sources into single strings
        collapse_uniprots_by_line(master_ortho_gaf)

    clean(master_ortho_gaf, header_lines = 0, dedup=False, sort=True, sort_cols = [2,4])

    # Add all orthologs to Xenbase GAF
    xen_w_ortho = f'{output_dir}/xenbase.plus.orthologs.gaf'
    combine_annotations(xenbase_gaf, master_ortho_gaf, xen_w_ortho)
    clean(xen_w_ortho, header_lines = 26, dedup=False, sort=True, sort_cols = [2,6,4])

    # Compress final files, remove uncompressed version
    gzip_file(xen_w_ortho, delete_input=False)

    print("Finished!\n")

# ------------------------------------- Supporting Functions --------------------------------------

# FUNCTION: Create uniprot/ncbi dictionaries to map xenopus annotation data to orthologs
# Uses flags to avoid populating maps when not required
def populate_maps(gaf_in, ortho_lookup=None, species=None, split=False):
    
    # Get NCBI ID from DB_Xref string in GPI file
    def extract_ncbi_id(dbxref):
        match = re.search(r'NCBI_Gene:(\d+)', dbxref)
        return match.group(1) if match else None

    # Add GO IDs to uniprot map:
    def add_go_ids(gaf, uniprot_map):
        with open(gaf, 'r', encoding=encoding) as gaf_in:
            reader = csv.reader(gaf_in, delimiter='\t')
            for fields in reader:
                if not fields or fields[0].startswith('!'):
                    continue

                uniprot_id = fields[1]
                go_id = fields[4]

                if uniprot_id not in uniprot_map:
                    continue

                entry = uniprot_map[uniprot_id]
                
                if isinstance(entry, str):
                    uniprot_map[uniprot_id] = {
                        "ncbi_id": entry,
                        "go_ids": [go_id]
                    }
                elif isinstance(entry, dict):
                    entry.setdefault("go_ids", [])
                    if go_id not in entry["go_ids"]:
                        entry["go_ids"].append(go_id)
                else:
                    print(f"Unexpected value type for uniprot key {uniprot_id}: {type(entry)}")

        # Remove all uniprots that dont have associated go ids (ie. not found in gaf file)
        keys_to_drop = [
            up_id for up_id, entry in uniprot_map.items()
            if not isinstance(entry, dict) or not entry.get("go_ids")
        ]
        for up_id in keys_to_drop:
            del uniprot_map[up_id]

    # Prints records in map (default 10)
    def print_map(dict_name, num_records=10):
        for i, (key, info) in enumerate(dict_name.items()):
            if i >= num_records:
                break
            print(f"{key}: {info}")

    # Prepare xenbase <-> uniprot/ncbi maps (used in both xenbase gaf creation & adding ortholog annotations)
    if not species and not split:
        gpi = os.path.join(xb_dir, 'Xenbase.gpi')

        # Populate xenopus uniprot and xenopus ncbi -> xenbase info maps
        with open(gpi, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            for fields in reader:
                if not fields or fields[0].startswith('!'):
                    continue
                
                xenbase = fields[1]
                DB_symbol = fields[2]
                DB_Object_Name = fields[3]
                DB_Object_Synonym = fields[4]
                object_type = fields[5]
                dbxref = fields[8]
                uniprots = re.findall(r'Uni[Pp]rotKB:(.+?)(?:\||$)', dbxref)
                ncbi_id = extract_ncbi_id(dbxref)
                
                for uniprot in uniprots:
                    xen_uniprot_map.setdefault(uniprot, {})
                    xen_uniprot_map[uniprot].update({
                        "ncbi_id": ncbi_id,
                        "xenbase_id": xenbase,
                        "symbol": DB_symbol,
                        "name": DB_Object_Name,
                        "synonyms": DB_Object_Synonym,
                        "object_type": object_type,
                    })

                xen_ncbi_map.setdefault(ncbi_id, {})
                xen_ncbi_map[ncbi_id].update({
                    "xenbase_id": xenbase,
                    "symbol": DB_symbol,
                    "name": DB_Object_Name,
                    "synonyms": DB_Object_Synonym,
                    "object_type": object_type,
                    "uniprots": uniprots
                })

            print(f"\nLoaded {len(xen_uniprot_map)} Xenopus UniProt IDs -> Xenbase info into map")
            #print_map(xen_uniprot_map)

            print(f"\nLoaded {len(xen_ncbi_map)} Xenopus NCBI IDs -> Xenbase info into map")
            #print_map(xen_ncbi_map)

        # Add GO IDs and drop uniprots not in xenopus gaf
        add_go_ids(gaf_in, xen_uniprot_map)
        #print_map(xen_uniprot_map)

    # Prepare xenbase <-> ortholog NCBI maps
    if ortho_lookup:
        xenopus_orthologs = os.path.join(ncbi_map_dir, 'Xenopus_NCBI_Orthologs.tsv')

        # Populate xenbase ortholog lookup map:
        with open(xenopus_orthologs, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            header = next(reader)
            species_keys = [key.capitalize() for key in header[2:8]]

            for fields in reader:
                xen_ncbi_id = fields[1]
                if not xen_ncbi_id:
                    continue

                species_map = dict(zip(species_keys, fields[2:8]))

                # Xenopus NCBI -> Ortho Species -> Ortho NCBI
                xen_to_ortho_lookup[xen_ncbi_id] = species_map

                # Ortho Species -> Ortho NCBI -> Xenopus NCBI
                for ortho_species, ortho_ncbi in species_map.items():
                    if not ortho_ncbi:
                        continue

                    if ortho_species not in ortho_to_xen_lookup:
                        ortho_to_xen_lookup[ortho_species] = {}

                    ortho_to_xen_lookup[ortho_species][ortho_ncbi] = xen_ncbi_id

            print(f"\nLoaded {len(xen_to_ortho_lookup)} Xenopus NCBI IDs -> Ortholog NCBI IDs into map")
            #print_map(xen_to_ortho_lookup)
            print(f"\nLoaded {len(ortho_to_xen_lookup['Human'])} Ortholog NCBI IDs (Ex. Human) -> Xenopus NCBI IDs into map")
            #print_map(ortho_to_xen_lookup['Human'])

    # Prepare ortholog uniprot map (NOTE: Dependent on ortho lookup maps; repopulated for each species)
    if species:
        ortho_mapping_file = os.path.join(ncbi_map_dir, f"{species}_NCBI_Mapping.tsv")

        if not os.path.exists(ortho_mapping_file):
            print(f"ERROR: No NCBI mapping file found for {species} at '{ortho_mapping_file}'")
            print(f"Skipping {species} Uniprot & NCBI map creation")
            return
        
        # Populate ortholog uniprot & ncbi maps:
        with open(ortho_mapping_file, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            orthologus_ncbi_ids = {
                entry[species] for entry in xen_to_ortho_lookup.values() if entry[species]
            }

            for fields in reader:
                uniprot_id = fields[0]
                ncbi_id = fields[2]

                if ncbi_id not in orthologus_ncbi_ids:
                    continue        # If ncbi id is not found in ortholog lookup table, there is no xenopus match

                ortho_uniprot_map[uniprot_id] = ncbi_id

            print(f"\nLoaded {len(ortho_uniprot_map)} {species} UniProt IDs -> NCBI IDs into map")
            #print_map(ortho_uniprot_map)

        # Add GO IDs and drop uniprots not in ortholog gaf (from uniprot map)
        add_go_ids(gaf_in, ortho_uniprot_map)
        #print_map(ortho_uniprot_map)

    # Prepare maps for splitting xenbase gaf by xenopus species
    if split:
        xb_genepage_to_geneid = os.path.join(xb_dir, 'Xenbase_Genepage_To_GeneId.txt')

        # Populate x.trop & x.laev gene id maps
        with open(xb_genepage_to_geneid, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            header = next(reader)
            for fields in reader:
                umbrella = fields[0]
                xtrop = fields[2]
                xlaevl = fields[4]
                xlaevs = fields[6]

                if xtrop:
                    xtrop_map[xtrop] = umbrella
                if xlaevl:
                    xlaev_map[xlaevl] = umbrella
                if xlaevs:
                    xlaev_map[xlaevs] = umbrella
            
            print(f"\nLoaded {len(xtrop_map)} X.trop IDs -> genepage IDs into map")
            #print_map(xtrop_map)
            print(f"Loaded {len(xlaev_map)} X.laev IDs -> genepage IDs into map")
            #print_map(xlaev_map)

# FUNCTION: Filter GAF and return path of new filtered GAF
def filter_gaf(gaf_in, create_copy=True, filter_values=None, filter_col=None):
    print(f"\nFiltering {os.path.basename(gaf_in)} on {gaf_columns[filter_col]}...")
    gaf_name = os.path.basename(gaf_in)[:-4]   # remove extension
    gaf_out = os.path.join(gaf_dir, gaf_name + "_filtered.tmp")

    with open(gaf_in, 'r', encoding=encoding) as f_in, open(gaf_out, 'w', encoding=encoding) as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        
        for fields in reader:
            if fields[0].startswith('!'):
                #print("Skipping header line...")
                writer.writerow(fields)
                continue

            if filter_col and filter_values:
                if fields[filter_col] not in filter_values:
                    continue

            writer.writerow(fields)

    if not create_copy:
        shutil.move(gaf_out, gaf_in)
        gaf_out = gaf_in

    return gaf_out      # returns path of filtered gaf for use later if needed

# FUNCTION: Write header for matched annotations file
def write_gaf_header(gaf_output):
    header = textwrap.dedent(f"""\
        !gaf-version: 2.2
        !Project_name: Xenbase 
        !URL: http://www.xenbase.org/
        !Contact Email: xenbase@cchmc.org
        !Funding: NICHD of US National Institutes of Health
        !date: {date.today()}
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
        !
    """)
    with open(gaf_output, 'w') as gaf:
        gaf.write(header)

# FUNCTION: Parse GAF and pass fields into matching function (Adapted from original GOA_parsing perl script)
def parse_gaf(gaf_in, process_func):
    filename = os.path.basename(gaf_in)[:-4]

    with open(gaf_in, 'r', encoding=encoding) as gaf:
        reader = csv.reader(gaf, delimiter='\t')
        for i, fields in enumerate(reader, start=1):
            if not fields or fields[0].startswith('!'):
                continue
            if len(fields) < 17:
                continue

            records = {
                "db": fields[0],
                "object_id": fields[1],
                "symbol": fields[2],
                "qualifier": fields[3],
                "go_id": fields[4],
                "db_ref": fields[5],
                "evidence": fields[6],
                "with_from": fields[7],
                "aspect": fields[8],
                "object_name": fields[9],
                "object_synonyms": fields[10],
                "object_type": fields[11],
                "taxon": fields[12],
                "date": fields[13],
                "assigned_by": fields[14],
                "annotation_extension": fields[15],
                "gene_product_id": fields[16]
            }

            process_func(records)

# FUNCTION: Create Xenbase matched & unmatched GAFs and provenance file ('matched' + original db and uniprot info) (Adapted from original GOA_parsing perl script)
def create_xenbase_gaf(gaf_in, output_dir, zip=True):
    def wrapper(fields):
        match_to_xenbase(fields, matched, provenance, unmatched)

    # FUNCTION: Find rows in Xenopus GAF where uniprot ID exists in xen_uniprot_map (Adapted from original GOA_parsing perl script)
    # Matched & provenance output will ONLY contain annotations where a xenbase gene ID was found
    # Unmatched output will ONLY contain annotations where NO xenbase genes were found
    def match_to_xenbase(fields, matched, provenance=None, unmatched=None):
        uniprot_id = fields["object_id"]
        if uniprot_id in xen_uniprot_map:
            xb_data = xen_uniprot_map[uniprot_id]

            # Creates main gaf file of xenopus go annotations linked with xenbase ids
            # Subs:
            #   - db -> 'Xenbase'
            #   - uniprot id -> xb gene symbol
            #   - symbol -> symbol from xb
            #   - object name -> name from xb
            #   - object synonyms -> synonyms from xb
            with open(matched, 'a') as m:
                m.write(f"Xenbase\t{xb_data['xenbase_id']}\t{xb_data['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{fields['with_from']}\t{fields['aspect']}\t{xb_data['name']}\t{xb_data['synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")

            if provenance:                
                with open(provenance, 'a') as p:
                    p.write(f"Xenbase\t{xb_data['xenbase_id']}\t{xb_data['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{fields['with_from']}\t{fields['aspect']}\t{xb_data['name']}\t{xb_data['synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['db']}:{fields['object_id']}\n")
        else:
            if unmatched:
                with open(unmatched, 'a') as u:
                    u.write(f"{fields['db']}\t{fields['object_id']}\t{fields['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{fields['with_from']}\t{fields['aspect']}\t{fields['object_name']}\t{fields['object_synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")    
   
    print(f"\n------------------------ Creating Xenbase GAF ------------------------\n")

    matched = os.path.join(output_dir, f'xenbase.EBI.only.gaf')
    unmatched = os.path.join(output_dir, f'unmatched_annotations.tsv')
    provenance = os.path.join(output_dir, 'annotation_provenance.tsv')

    # Remove output files if they already exist (for clean write)
    for file in [matched, unmatched, provenance]:
        if file and os.path.exists(file):
            os.remove(file)

    # Write Xenbase header to output GAF
    write_gaf_header(matched)
    
    # Parse & match GAF with Xenbase GPI
    parse_gaf(gaf_in, wrapper)
    
    # Deduplicate on full line & sort by symbol then go id
    clean(matched, header_lines=26, dedup=True, sort=True, sort_cols=[2,4])                 # dedup_cols=[1,2,3,4,5,6,7]
    clean(provenance, dedup=True, sort=True, sort_cols=[2,4])

    # Optional: zip output
    if zip:
        gzip_file(matched)
        gzip_file(provenance, delete_input=True)
        gzip_file(unmatched, delete_input=True)

# FUNCTION: Split xenbase GAF into x.trop & x.laev specific files
def split_xenbase_gaf(xenbase_gaf, output_dir, zip=True):
    def wrapper(fields):
        split_by_xenopus(fields, xtrop_only, xlaev_only)

    def split_by_xenopus(fields, xtrop_file, xlaev_file):
        xb_gene_id = fields["object_id"]

        ordered_fields = [
        "db", "object_id", "symbol", "qualifier", "go_id", "db_ref", "evidence",
        "with_from", "aspect", "object_name", "object_synonyms", "object_type",
        "taxon", "date", "assigned_by", "annotation_extension", "gene_product_id"
        ]

        line = "\t".join(fields[key] for key in ordered_fields)

        if xb_gene_id in xtrop_map:
            with open(xtrop_file, 'a') as xtrop:
                xtrop.write(f"{line}\n")
        elif xb_gene_id in xlaev_map:
            with open(xlaev_file, 'a') as xlaev:
                xlaev.write(f"{line}\n")

    xtrop_only = os.path.join(output_dir, f'xenbase.xtrop.only.gaf')
    xlaev_only = os.path.join(output_dir, f'xenbase.xlaev.only.gaf')

    # Remove output files if they already exist
    for file in [xtrop_only, xlaev_only]:
        if file and os.path.exists(file):
            os.remove(file)

    # Write Xenbase header to output GAFs
    write_gaf_header(xtrop_only)
    write_gaf_header(xlaev_only)
    
    # Parse & match GAF with Xenbase GPI
    parse_gaf(xenbase_gaf, wrapper)

    # Optional: zip output
    if zip:
        gzip_file(xtrop_only, delete_input=True)
        gzip_file(xlaev_only, delete_input=True)

# FUNCTION: Sub xenbase info into matching ortholog annotations 
def modify_ortho_annotations(ortho_gaf, output_dir, species):
    def wrapper(fields):
        match_ortho_to_xenbase(fields, matched, species)

    # FUNCTION: Find rows in ortholog GAF where uniprot can be mapped to a xenbase ID. Sub in xenbase info and change annotation evidence to "ISO"
    def match_ortho_to_xenbase(fields, matched, species):
        xen_ncbi = None
        uniprot_id = fields["object_id"]

        if uniprot_id in ortho_uniprot_map:
            ortho_ncbi = ortho_uniprot_map[uniprot_id]["ncbi_id"]
            if ortho_ncbi in ortho_to_xen_lookup[species]:
                xen_ncbi = ortho_to_xen_lookup[species][ortho_ncbi]

        if xen_ncbi and xen_ncbi in xen_ncbi_map:
            PLACEHOLDER_REF = "GO_REF:PLACEHOLDER"
            TAXON_ID = "taxon:8364"
            EVIDENCE = "ISO"
            #SOURCE = species
            SOURCE = "Xenbase"
            DATE = date.today().strftime("%Y%m%d")
            xb_data = xen_ncbi_map[xen_ncbi]

            with open(matched, 'a') as m:
                m.write(f"Xenbase\t{xb_data['xenbase_id']}\t{xb_data['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{PLACEHOLDER_REF}\t{EVIDENCE}\tUniProtKB:{uniprot_id}\t{fields['aspect']}\t{xb_data['name']}\t{xb_data['synonyms']}\t{xb_data['object_type']}\t{TAXON_ID}\t{DATE}\t{SOURCE}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")

    matched = os.path.join(output_dir, f'ortho-gafs/Xenbase_from_{species}.gaf')

    # Remove matched file if it already exists (for clean write)
    if os.path.exists(matched):
        os.remove(matched)

    print(f"\nParsing {species} GOA for Xenbase matches...\n")
    parse_gaf(ortho_gaf, wrapper)

    # Remove duplicate lines and collapse uniprots for matching annotations into 1 line
    clean(matched, header_lines=0, dedup=True, dedup_cols=[1,4,7])
    collapse_uniprots_by_line(matched)
    
    # Return file path of modified ortho gaf
    return matched

# FUNCTION: Collapses annotation lines that only differ on the uniprot/source. Uniprots and sources are concatenated into a string
def collapse_uniprots_by_line(input_gaf):
    grouped = {}
    basename = os.path.basename(input_gaf)
    temp_file = f'{input_gaf}.tmp'

    with open(input_gaf, 'r', encoding=encoding) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            cols = line.strip().split('\t')
            
            # Build key, excluding With/From (7) and Source (14) columns
            key = tuple(cols[:7] + cols[8:14] + cols[15:])

            if key not in grouped:
                grouped[key] = {
                    'uniprots': [],
                    'sources': [],
                    'template': cols
                }

            if cols[7] not in grouped[key]['uniprots']:
                grouped[key]['uniprots'].append(cols[7])

            if cols[14] not in grouped[key]['sources']:
                grouped[key]['sources'].append(cols[14])

    with open(temp_file, 'w', encoding=encoding) as temp:
        for group in grouped.values():
            new_row = group['template'][:]
            new_row[7] = '|'.join(group['uniprots'])
            new_row[14] = '/'.join(group['sources'])
            temp.write('\t'.join(new_row) + '\n')

    shutil.move(temp_file, input_gaf)
    print(f"Concatenated uniprots in {basename}")
    #clean(input_gaf, header_lines = 0, dedup=False, sort=True, sort_cols = [2,4])

# FUNCTION: Adds annotations between 2 (GAF) files together
def combine_annotations(file_1, file_2, output_file):
    with open(file_1, 'r', encoding=encoding) as f1:
        file_1_lines = f1.readlines()

    with open(file_2, 'r', encoding=encoding) as f2:
        file_2_lines = f2.readlines()

    with open(output_file, 'w', encoding=encoding) as output:
        output.writelines(file_1_lines)
        output.writelines(file_2_lines)

# FUNCTION: Clean output files by deduplicating and/or sorting
def clean(filepath, header_lines = 0, dedup=True, dedup_cols=None, sort=False, sort_cols = None):
    filename = os.path.basename(filepath)

    # FUNCTION: removes duplicate lines either based on a full line match or a partial line match in a given set of columns
    # Default deduplication: full line match; for partial line match, set dedup_cols = [col_a, col_b, etc.] (ex. dedup_cols = [1,4] will remove duplicates based on gene id + go term pairs only)
    def deduplicate(lines, dedup_cols):
        seen = set()
        unique = []
        removed_count = 0

        for line in lines:
            if line.startswith('!'):
                unique.append(line)
                continue

            cols = line.split('\t')
            key = tuple(cols[i] for i in dedup_cols if i < len(cols)) if dedup_cols else line

            if key not in seen:
                unique.append(line)
                seen.add(key)
            else:
                removed_count += 1

        print(f"{removed_count} duplicate lines removed")
        return unique

    # FUNCTION: sorts lines in file by a given set of column(s,) giving priority to columns in the order they appear in the list (left to right)
    # Default sorting: file is not sorted; for sorting on multiple cols, set sort_cols = [col_a, col_b, etc.] (ex. sort_cols = [2, 4] will sort gaf by gene symbols, then go ids)    
    def sort_lines(lines, sort_cols):
        def sort_key(line):
            cols = line.split('\t')
            return tuple(cols[i].lower() for i in sort_cols if i < len(cols))
        return sorted(lines, key=sort_key)
        
    with open(filepath, "r") as f:
        lines = f.readlines()

    if header_lines != 0:
        # Save header
        header = lines[:header_lines] if isinstance(header_lines, int) else header_lines(lines) #26 for gaf, 1 for single header
        data_lines = lines[header_lines:]
    else:
        header = []
        data_lines = lines

    if dedup:
        if dedup_cols != None:
            dedup_col_names = [gaf_columns[i] for i in dedup_cols]
            print(f"Deduplicating {filename} by column(s): {", ".join(dedup_col_names)} ...")
        else:
            print("Deduplicating " + filename + " by full line...")
        data_lines = deduplicate(data_lines, dedup_cols)

    if sort:
        sort_col_names = [gaf_columns[i] for i in sort_cols]
        print(f"Sorting {filename} on column(s): {", ".join(sort_col_names)} ...")
        data_lines = sort_lines(data_lines, sort_cols)

    with open(filepath, "w") as f:
        f.writelines(header)
        f.writelines(data_lines)

    #print(f"{filename} Cleaned!\n")
    print("")

# FUNCTION: Remove temp files from given folder
def remove_tmp_files(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith(".tmp"):
            file_path = os.path.join(folder_path, filename)
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

# FUNCTION: Compress files into .gz format
def gzip_file(input_path, delete_input=False):
    output_path = f'{input_path}.gz'
    with open(input_path, 'rb') as f_in:
        with gzip.open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if delete_input:
        os.remove(input_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--create_gaf', action='store_true')
    parser.add_argument('--add_orthos', nargs='?', help='List of ortholog species to process')
    parser.add_argument('--date', nargs='?', help='Date of download for files used')
    parser.add_argument('--log', action='store_true')
    args = parser.parse_args()

    # ---------------------------- Variable Setup ----------------------------
    HOME = os.path.expanduser("~/xenbase-gaf-pipeline")

    # Define folder paths
    global input_dir, gaf_dir, ncbi_map_dir, xb_dir, output_dir
    input_dir = os.path.join(HOME, "input-files")
    gaf_dir = os.path.join(input_dir, "goa-gafs")
    ncbi_map_dir = os.path.join(input_dir, "ncbi-maps")
    xb_dir = os.path.join(input_dir, "xenbase-files")
    output_dir = os.path.join(HOME, "output-files")

    # Ensure folders exist
    for folder in [input_dir, gaf_dir, ncbi_map_dir, xb_dir, output_dir]:
        os.makedirs(folder, exist_ok=True)

    # Redirect output to log file
    if args.log:
        log_path = os.path.join(HOME, "output-files/script-logs/goa_parsing.log")
        if args.create_gaf:
            log = open(log_path, "wt")
        elif args.add_orthos:
            log = open(log_path, "at")
        sys.stdout = log
        sys.sterr = log
    
    # Set gaf properties
    global encoding, gaf_columns, allowed_codes
    encoding = "utf-8"
    gaf_columns = ["DB", "DB Object ID", "DB Object Symbol", "Qualifier", "GO ID",
      "DB:Reference", "Evidence Code", "With/From", "Aspect", "DB Object Name",
      "DB Object Synonym", "DB Object Type", "Taxon", "Date", "Assigned By",
      "Annotation Extension", "Gene Product ID"]
    allowed_codes = {"EXP", "IDA", "IMP", "IGI", "IEP", "ISS", "TAS", "IC"}

    # Define maps
    global xen_uniprot_map, xen_ncbi_map, xtrop_map, xlaev_map
    global xen_to_ortho_lookup, ortho_to_xen_lookup, ortho_uniprot_map
    xen_uniprot_map = {}                # key = xenopus uniprot id: contains ncbi id, xenbase gene id, symbol, name, synonyms, & go ids
    xen_ncbi_map = {}                   # key = xenopus ncbi id: contains xenbase gene id, symbol, name, synonyms & uniprots
    xtrop_map = {}                      # key = x.trop gene id; contains xenbase genepage id
    xlaev_map = {}                      # key = x.laev(.L/.S) gene id; contains xenbase genepage id
    xen_to_ortho_lookup = {}            # key = ortho species, ortholog ncbi id: contains xenopus ncbi ids
    ortho_to_xen_lookup = {}            # key = xenopus ncbi id, ortho species: contains ortholog ncbi ids
    ortho_uniprot_map = {}              # key = ortholog uniprot id: contains ncbi id & go ids

    # Set file download date to today if no date flag was provided
    dl_date = args.date if args.date else datetime.today().strftime('%Y-%m-%d')
    print(f"\nDate & time of script execution: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")
    
    if args.create_gaf:
        main_gaf(dl_date)
    elif args.add_orthos:
        ortho_species = ast.literal_eval(args.add_orthos)
        main_ortho(dl_date, ortho_species)
    else:
        print("""Error: Must specify processing option:
        --create_gaf    -> Use to create xenbase GAF from EBI's GOA files
        --add_orthos    -> Use to add ortholog annotations to xenbase gaf (NOTE: use after noctua pipeline)""")