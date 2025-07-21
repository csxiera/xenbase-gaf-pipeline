import os
import gzip
import shutil
import requests
import re
import csv
import textwrap
from datetime import date
from io import StringIO

# NOTE: Use get_files.sh to download required files if needed
# Original GOA_parsing perl script adapted into create_xenbase_gaf, parse_gaf, and match_to_xenbase functions

# *** Need to add GOLR annotations as well? - if gaf format is the same as GOA gaf, should be able to use same functions for processing

def main():
    # ---------------------------- Setup ----------------------------
    HOME = os.path.expanduser("~/xenbase-gaf-pipeline")
    #script_dir = os.path.dirname(os.path.abspath(__file__))

    # folder paths
    global input_dir, gaf_dir, ncbi_map_dir, output_dir
    input_dir = os.path.join(HOME, "input-files")
    gaf_dir = os.path.join(input_dir, "gafs")
    ncbi_map_dir = os.path.join(input_dir, "ncbi-map")
    output_dir = os.path.join(HOME, "output-files")
    
    global encoding
    encoding = "latin1" # To handle accented characters

    global gaf_columns, allowed_codes
    gaf_columns = ["DB", "DB Object ID", "DB Object Symbol", "Qualifier", "GO ID",
      "DB:Reference", "Evidence Code", "With/From", "Aspect", "DB Object Name",
      "DB Object Synonym", "DB Object Type", "Taxon", "Date", "Assigned By",
      "Annotation Extension", "Gene Product ID"]
    allowed_codes = {"EXP", "IDA", "IMP", "IGI", "IEP", "ISS", "TAS", "IC"}

    # Maps
    global xen_uniprot_map, xen_ncbi_map, ortholog_lookup, ortho_uniprot_map, ortho_ncbi_to_uniprot
    xen_uniprot_map = {}                # key = xenopus uniprot id: contains ncbi id, xenbase id, symbol, name, synonyms, & go ids
    xen_ncbi_map = {}                   # key = xenopus ncbi ids: contains xenbase id, symbol, name, synonyms, etc.
    ortholog_lookup = {}                # key = ortholog ncbi ids: contains xenopus ncbi ids
    ortho_uniprot_map = {}              # key = ortholog uniprot id: contains ncbi id & go ids
    ortho_ncbi_to_uniprot = {}          # key = ortholog ncbi id: contains uniprot id

    ortho_species_list = ["Human", "Mouse", "Rat", "Chicken", "Zebrafish", "Drosophila"]
    # FIX!!: want to pass this info in maybe??
    xenopus_gaf = os.path.join(gaf_dir, 'Xenopus.GOA.Extracted.2025-07-18.gaf')
    ortho_gaf = os.path.join(gaf_dir, 'Human.GOA.Extracted.2025-07-18.gaf')

    '''xenopus_gaf = os.path.join(gaf_dir, 'Xenopus.GOA.TEST.gaf')
    ortho_gaf = os.path.join(gaf_dir, 'Human.GOA.TEST.gaf')'''
    ortho_ncbi_mapping = os.path.join(ncbi_map_dir, 'Human_NCBI_Mapping.tsv')

    # ---------------------------- Data Processing ----------------------------

    sub_ortholog_uniprots = False    # Set to false if Xenbase GAF without ortholog uniprots is needed

    if sub_ortholog_uniprots:
        # Filter ortholog GAF file to remove invalid evidence codes
        filtered_ortho_gaf = filter_gaf(ortho_gaf, create_copy=True, filter_col=6, filter_values=allowed_codes)

        # Populates maps needed to sub ortholog uniprot IDs into Xenbase GAF with/from col
        populate_maps(xenopus_gaf, 'Xenopus')
        populate_maps(filtered_ortho_gaf, 'Human', ortho_ncbi_mapping)

        # Filter ortholog gaf to only include annotations with xenopus ncbi gene id matches
        filtered_ortho_gaf = filter_gaf(filtered_ortho_gaf, create_copy=False, filter_col=1, filter_values=ortho_uniprot_map.keys())
        clean(filtered_ortho_gaf, sort=True, sort_cols=[2,1,4]) # sort by symbol, then uniprot, then go id

        # Replace db, object_id, etc in human gaf with mapped xb info

        # Create Xenbase GAF from GPI & add in human annotations
        create_xenbase_gaf(xenopus_gaf, output_dir)
    
    else:
        # Modify GAF (from EBI) with xenbase gene information ONLY (no ortholog uniprots in with/from col)
        populate_maps(xenopus_gaf, 'Xenopus')
        create_xenbase_gaf(xenopus_gaf, output_dir)

    print("\nFinished.\n")

# ------------------------------- Functions -------------------------------
def set_files(ortho_species_list):

    pass

# FUNCTION: Create uniprot/symbol dictionary to map xenopus annotation data to orthologs
def populate_maps(gaf_in, species, ortho_mapping_file=None):
    
    # FUNCTION: Get NCBI ID from DB_Xref string in GPI file
    def extract_ncbi_id(dbxref):
        match = re.search(r'NCBI_Gene:(\d+)', dbxref)
        return match.group(1) if match else None

    # FUNCTION: Add GO IDs to uniprot map:
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

        # Keep record of uniprot in map but initialize blank go id list
        """for up_id, entry in uniprot_map.items():
            if isinstance(entry, str):
                uniprot_map[up_id] = {
                    "ncbi_id": entry,
                    "go_ids": []
                }
            elif isinstance(entry, dict) and "go_ids" not in entry:
                entry.setdefault("go_ids", [])"""

        print(f"\nAdded GO IDs to {species} Uniprot map")

    # FUNCTION: Prints records in map (default 10)
    def print_map(dict_name, num_records=10):
        for i, (key, info) in enumerate(dict_name.items()):
            if i >= num_records:
                break
            print(f"{key}: {info}")

    # Prepare xenopus map
    if species == 'Xenopus': 
        gpi = os.path.join(input_dir, 'Xenbase.gpi')

        # Populate xenopus uniprot -> xenbase info map
        # FIX!!: Create xenopus ncbi id -> xenbase info map
        with open(gpi, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            for fields in reader:
                if not fields or fields[0].startswith('!'):
                    continue
                
                xenbase = fields[1]
                DB_symbol = fields[2]
                DB_Object_Name = fields[3]
                DB_Object_Synonym = fields[4]
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
                        "synonyms": DB_Object_Synonym
                    })

            print(f"\nLoaded {len(xen_uniprot_map)} Xenopus UniProt IDs -> Xenbase info into map")
            #print_map(xen_uniprot_map)

        # Add GO IDs and drop any uniprots without them (ie. not found in xenopus gaf)
        add_go_ids(gaf_in, xen_uniprot_map)
        #print_map(xen_uniprot_map)

    # Prepare ortholog maps (if ortholog file/species provided)
    elif species != 'Xenopus' and ortho_mapping_file:
        xenopus_orthologs = os.path.join(ncbi_map_dir, 'Xenopus_NCBI_Orthologs.tsv')

        # Populate ortholog lookup map:
        # FIX!!: change ortholog lookup map to be [species][ortho_ncbi] = xen_ncbi
        with open(xenopus_orthologs, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            header = next(reader)
            species_keys = [key.capitalize() for key in header[2:8]]
            print(species_keys)

            for fields in reader:
                xen_ncbi_id = fields[1]
                if not xen_ncbi_id:
                    continue

                species_map = dict(zip(species_keys, fields[2:8]))
                ortholog_lookup[xen_ncbi_id] = species_map

            print(f"\nLoaded {len(ortholog_lookup)} Xenopus NCBI gene IDs -> orthologs into map")
            #print_map(ortholog_lookup)

        # Populate ortholog uniprot & ncbi maps:
        with open(ortho_mapping_file, 'r', encoding=encoding) as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            orthologus_ncbi_ids = {
                entry[species] for entry in ortholog_lookup.values() if entry[species]
            }

            for fields in reader:
                uniprot_id = fields[0]
                ncbi_id = fields[2]

                if ncbi_id not in orthologus_ncbi_ids:
                    continue        # If human ncbi id is not found in ortholog lookup table, there is no xenopus match

                ortho_uniprot_map[uniprot_id] = ncbi_id

                if ncbi_id in ortho_ncbi_to_uniprot:
                    ortho_ncbi_to_uniprot[ncbi_id].append(uniprot_id)
                else:
                    ortho_ncbi_to_uniprot[ncbi_id] = [uniprot_id]

            print(f"Loaded {len(ortho_uniprot_map)} {species} UniProt IDs -> NCBI gene IDs into map")
            #print_map(ortho_uniprot_map)
            print(f"Loaded {len(ortho_ncbi_to_uniprot)} {species} NCBI gene IDs -> Uniprot IDs into map")
            #print_map(ortho_ncbi_to_uniprot)

        # Add GO IDs and drop any uniprots without them (ie. those not found in ortholog gaf)
        add_go_ids(gaf_in, ortho_uniprot_map)
        #print_map(ortho_uniprot_map)

        # Filter out uniprots that arent found in ortholog gaf
        for ncbi_id in list(ortho_ncbi_to_uniprot.keys()):
            uniprot_list = ortho_ncbi_to_uniprot[ncbi_id]
            valid_uniprots = []

            for uniprot in uniprot_list:
                if uniprot in ortho_uniprot_map.keys():
                    valid_uniprots.append(uniprot)
            
            if valid_uniprots:
                ortho_ncbi_to_uniprot[ncbi_id] = valid_uniprots
            else:
                del ortho_ncbi_to_uniprot[ncbi_id]
        #print_map(ortho_ncbi_to_uniprot)

# FUNCTION: Filter gaf file and return path of new filtered gaf
def filter_gaf(gaf_in, create_copy=True, filter_values=None, filter_col=None):
    print(f"Filtering {os.path.basename(gaf_in)} on {gaf_columns[filter_col]} ...")
    gaf_name = os.path.basename(gaf_in)[:-4]   # remove extension
    gaf_out = os.path.join(gaf_dir, gaf_name + "_filtered.gaf")

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

    print("Filtered output saved to " + gaf_out + "\n")
    return gaf_out      # returns path of filtered gaf for use later if needed

# FUNCTION: Sub xenbase info into matching ortholog annotations 
def modify_ortho_annotations(ortho_gaf, output_dir, species):
    def wrapper(fields):
        match_to_xenbase(fields, matched, species)

    matched = os.path.join(output_dir, f'ortho-gafs/Xenbase_from_{species}.gaf')
    parse_gaf(ortho_gaf, wrapper)

# FUNCTION: Create Xenbase matched & unmatched GAFs and provenance file ('matched' + original db and uniprot info) (Adapted from original GOA_parsing perl script)
def create_xenbase_gaf(gaf_in, output_dir):
    def wrapper(fields):
        match_to_xenbase(fields, matched, provenance, unmatched)

    filename = os.path.basename(gaf_in)[:-4]
    
    matched = os.path.join(output_dir, f'Xenbase.gaf')
    unmatched = os.path.join(output_dir, f'unmatched_annotations.{filename}.gaf')
    provenance = os.path.join(output_dir, 'annotation_provenance.tmp')

    # Remove output files if they already exist (for clean write)
    for file in [matched, unmatched, provenance]:
        if file and os.path.exists(file):
            os.remove(file)

    # Write Xenbase header to output GAF
    write_gaf_header(matched)
    
    # Parse & match GAF with Xenbase GPI
    parse_gaf(gaf_in, wrapper)
    
    # Deduplicate on object id, qualifier, go id, evidence code, and with/from columns
    # Sort by symbol then go id
    clean(matched, header_lines=26, dedup=True, dedup_cols=[1,3,4,6,7], sort=True, sort_cols=[2,4])
    clean(provenance, dedup=True, dedup_cols=[1,3,4,6,7], sort=True, sort_cols=[2,4])

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
    print(f"\nFinding Xenbase matches in {filename}...\n")

    with open(gaf_in, 'r', encoding=encoding) as gaf:
        reader = csv.reader(gaf, delimiter='\t')
        for i, fields in enumerate(reader, start=1):
            if not fields or fields[0].startswith('!'):
                continue
            if len(fields) < 17:
                continue

            records = {
                "db": fields[0],
                "uniprot_id": fields[1],
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

# FUNCTION: Find rows in GAF where uniprot ID exists in xen_uniprot_map (Adapted from original GOA_parsing perl script)
# Matched & provenance output will ONLY contain annotations where a xenbase gene ID was found
# Unmatched output will ONLY contain annotations where NO xenbase genes were found
def match_to_xenbase(fields, matched, provenance=None, unmatched=None):
    #print(f"Uniprot passed through fields: {fields["uniprot_id"]}")
    uniprot_id = fields["uniprot_id"]
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
                p.write(f"Xenbase\t{xb_data['xenbase_id']}\t{xb_data['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{fields['with_from']}\t{fields['aspect']}\t{xb_data['name']}\t{xb_data['synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['db']}:{fields['uniprot_id']}\n")
    else:
        if unmatched:
            with open(unmatched, 'a') as u:
                u.write(f"{fields['db']}:{fields['uniprot_id']}\t{fields['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{fields['with_from']}\t{fields['aspect']}\t{fields['object_name']}\t{fields['object_synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")    

# FUNCTION: Sub xenbase info into ortholog gaf annotations
def match_ortho_to_xenbase(fields, matched, species):
    uniprot_id = fields["uniprot_id"]

    if uniprot_id in ortho_uniprot_map:
        ncbi_id = ortho_uniprot_map[uniprot_id]
        for xenopus_ncbi, species_map in ortholog_lookup.items():   #FIX!!: change ortholog lookup map to be [species][ortho_ncbi] = xen_ncbi
            human_ids = species_map.get("human", [])
            if ncbi_id in human_ids:
                pass
                # Use xenopus ncbi id to access xb info?

# FUNCTION: Add ortholog annotations (as ISO evidence types) to xenbase GAF (as new GAF)
def add_ortho_annotations(xenbase_gaf, output_dir, ortho_species):
    for species in ortho_species:
        ortho_gaf = os.path.join(output_dir, f"ortho-gafs/Xenbase_from_{species}.gaf")
        output_gaf = os.path.join(output_dir, f"Xenbase_w_{species}.gaf")

        with open(xenbase_gaf, 'r', encoding=encoding) as xen_file:
            xen_lines = xen_file.readlines()

        with open(ortho_gaf, 'r', encoding=encoding) as ortho_file:
            ortho_lines = ortho_file.readlines()

        with open(output_gaf, 'w', encoding=encoding) as output:
            output.writelines(xen_lines)
            output.writelines(ortho_lines)

# FUNCTION: Clean output files by deduplicating and/or sorting
# Default deduplication: file is deduplicated on the full line, but can be deduplicated on specific column(s) (ex. dedup_cols = [1,4] will remove duplicates based on gene id + go term pairs only)
# Default sorting: file is not sorted, but can be sorted by multiple columns (ex. sort_cols = [2] will sort gaf by gene symbols)
def clean(filepath, header_lines = 0, dedup=True, dedup_cols=None, sort=False, sort_cols = None):
    filename = os.path.basename(filepath)

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
            
    def sort(lines, sort_cols):
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
            print("Deduplicating " + filename + " by column(s): " + str(dedup_cols) + "...")
        else:
            print("Deduplicating " + filename + " by full line...")
        data_lines = deduplicate(data_lines, dedup_cols)

    if sort:
        print("Sorting " + filename + " on column(s): " + str(sort_cols))
        data_lines = sort(data_lines, sort_cols)

    with open(filepath, "w") as f:
        f.writelines(header)
        f.writelines(data_lines)


main()


# Logic for subbing ortho uniport into EXISTING xenopus gaf annotations (WRONG APPROACH)
'''
# FUNCTION: Modify Xenbase GAF to include ortholog uniports in with/from (when present), and create new Xenbase GAF that only includes annotations with orthologs
def create_xenbase_w_ortho_gaf(gaf_in, output_dir, ortho_species):
    def ortho_wrapper(fields):
        match_to_ortho(fields, unfiltered=xen_w_ortho_temp, filtered=ortho_filtered_temp, ortho_species=ortho_species, check_go_id=match_go_ids)

    def xenbase_wrapper_unfiltered(fields):
        match_to_xenbase(fields, matched=xenbase_w_ortho)

    def xenbase_wrapper_filtered(fields):
        match_to_xenbase(fields, matched=xenbase_ortho_filtered)

    xen_w_ortho_temp = os.path.join(output_dir, f'Xenopus_w_{ortho_species}.tmp')
    ortho_filtered_temp = os.path.join(output_dir, f'Xenopus_{ortho_species}_Filtered.tmp')
    xenbase_w_ortho = os.path.join(output_dir, f'Xenbase_w_{ortho_species}.gaf')
    xenbase_ortho_filtered = os.path.join(output_dir, f'Xenbase_{ortho_species}_Filtered.gaf')
    match_go_ids = False    # Set to True to ensure ortholog uniprot match is associated with same GO ID as xenopus annotation

    # Remove intermediate/final output files if they already exist (for clean write)
    for file in [xen_w_ortho_temp, ortho_filtered_temp, xenbase_w_ortho, xenbase_ortho_filtered]:
        if file and os.path.exists(file):
            os.remove(file)

    # Parse & match original xenopus GAF (from EBI) with orthologs
    parse_gaf(gaf_in, ortho_wrapper)

    # Write Xenbase header to final output GAF
    write_gaf_header(xenbase_w_ortho)

    # Parse and match xenopus+ortholog GAF with Xenbase GPI
    parse_gaf(xen_w_ortho_temp, xenbase_wrapper_unfiltered)
    parse_gaf(ortho_filtered_temp, xenbase_wrapper_filtered)

    # Deduplicate on object id, qualifier, go id, evidence code, and with/from columns
    # Sort by symbol then go id
    clean(xenbase_w_ortho, header_lines=26, dedup=True, dedup_cols=[1,3,4,6,7], sort=True, sort_cols=[2,4])
    clean(xenbase_ortho_filtered, dedup=True, dedup_cols=[1,3,4,6,7], sort=True, sort_cols=[2,4])

    # Remove intermediate files
    for file in [xen_w_ortho_temp, ortho_filtered_temp]:
        if file and os.path.exists(file):
            os.remove(file)

# FUNCTION: Find rows in GAF where Xenbase uniprot ID can be mapped to an ortholog uniprot ID
# Unfiltered output will contain ALL annotations from original gaf, only subbing with/from IF ortholog uniprot found
# Filtered output will ONLY contain annotations where an ortholog uniprot was found
# Ortholog GAF info stored in ortho_uniprot_map and ortho_ncbi_to_uniprot map
def match_to_ortho(fields, unfiltered, filtered, ortho_species, check_go_id=False):
    uniprot_id = fields["uniprot_id"]
    go_id = fields["go_id"]
    valid_uniprots = set()

    if uniprot_id in xen_uniprot_map:
        ncbi_id = xen_uniprot_map[uniprot_id]["ncbi_id"]
        if ncbi_id in ortholog_lookup:
            ortho_ncbi = ortholog_lookup[ncbi_id][ortho_species]
            if ortho_ncbi and ortho_ncbi in ortho_ncbi_to_uniprot:
                ortho_uniprots = ortho_ncbi_to_uniprot[ortho_ncbi]
                for ortho_uniprot in ortho_uniprots:
                    if not check_go_id:
                        valid_uniprots.add(f"UniprotKB:{ortho_uniprot}")
                    else:
                        ortho_go_ids = ortho_uniprot_map.get(ortho_uniprot, set())
                        if go_id in ortho_go_ids:
                            valid_uniprots.add(f"UniprotKB:{ortho_uniprot}")

    if valid_uniprots:
        with_from = "|".join(sorted(valid_uniprots))
        with open(filtered, 'a') as f:
                f.write(f"{fields['db']}\t{fields['uniprot_id']}\t{fields['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{with_from}\t{fields['aspect']}\t{fields['object_name']}\t{fields['object_synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")
        with open(unfiltered, 'a') as u:
                u.write(f"{fields['db']}\t{fields['uniprot_id']}\t{fields['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t{with_from}\t{fields['aspect']}\t{fields['object_name']}\t{fields['object_synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")

    else:
        with open(unfiltered, 'a') as u:
                u.write(f"{fields['db']}\t{fields['uniprot_id']}\t{fields['symbol']}\t{fields['qualifier']}\t{fields['go_id']}\t{fields['db_ref']}\t{fields['evidence']}\t \t{fields['aspect']}\t{fields['object_name']}\t{fields['object_synonyms']}\t{fields['object_type']}\t{fields['taxon']}\t{fields['date']}\t{fields['assigned_by']}\t{fields['annotation_extension']}\t{fields['gene_product_id']}\n")
'''





