import argparse
import sys
import requests
import os
import subprocess
import gzip
import shutil
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

# AUTHOR: C. Lenz
#
# SCRIPT FUNCTION: 
#   1. Download GOA GAF file and extract into species-specific gafs
#   2. Download NCBI -> uniprot maps for all potential ortholog species
#   3. Download xenbase files, such as GPI and NCBI ortholog mapping file
#
# NOTE: URL strings may change in the future! 
# Update EBI/uniprot/xenbase URL string in appropriate function if/when this happens
#
# Expected download times (from July 2025):
#   - goa_uniprot_all.gaf.gz (~4hrs, 21GB)
#   - species-specific gafs (~18 min each)
#   - remaining files (>1 min each)

# FUNCTION: Setup folders if they dont yet exist
def set_folders():
    paths = [
        input_dir,
        os.path.join(out_dir, "script-logs"),
        gaf_dir,
        os.path.join(gaf_dir, "original-goa"),
        ncbi_map_dir,
        xb_dir,
        out_dir
    ]
    for path in paths:
        os.makedirs(path, exist_ok=True)

# FUNCTION: Download GOA file with data for all species
# Contains EBI URL
def download_goa(args):
    if args.xen_download:
        xtrop_url = "https://ftp.ebi.ac.uk/pub/contrib/goa/gp_association.8364_Xenopus_tropicalis.gaf.gz"
        xtrop_path = os.path.join(gaf_dir, "original-goa", "8364_Xenopus_tropicalis.gaf.gz")
        download_w_progress(xtrop_url, output_path=xtrop_path)

        xlaev_url = "https://ftp.ebi.ac.uk/pub/contrib/goa/gp_association.8355_Xenopus_laevis.gaf.gz"
        xlaev_path = os.path.join(gaf_dir, "original-goa", "8355_Xenopus_laevis.gaf.gz")
        download_w_progress(xlaev_url, output_path=xlaev_path)

        xenopus_combined = os.path.join(gaf_dir, f"Xenopus.GOA.Curated.gaf")
        combine_annotations(xtrop_path, xlaev_path, xenopus_combined)

    elif args.ortho_download:
        print("Downloading GOA file from EBI...")
        goa_url = f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
        output_path = os.path.join(gaf_dir, "original-goa", GOA_FILENAME)
        # chunk size = <BYTE SIZE> (To modify number of bytes streamed in a chunk during download; default = 10*1024*1024 (10MB))
        download_w_progress(goa_url, output_path=output_path)

# FUNCTION: Extract GAF data for each species into seperate files
# NOTE: Extraction progress is not logged to get_files.log -> will output to terminal
def extract_from_goa(species):
    input_path = os.path.join(gaf_dir, "original-goa", GOA_FILENAME)
    if not os.path.exists(input_path):
        print(f"ERROR: '{input_path}' does not exist")
        return
    
    targets = species_taxon_map.keys() if species == "all" else [species.lower()]
    for target in targets:
        if target not in species_taxon_map:
            print(f"ERROR: '{target}' is not a valid species in taxon map")
            continue

        print(f"Extracting {target.capitalize()} annotations from '{GOA_FILENAME}'...")

        # Get taxon pattern to look for in GOA file (if multiple taxons for given target, join by "\|")
        taxon_list = species_taxon_map[target]
        taxon_pattern = r'\|'.join([f"taxon:{taxon}" for taxon in taxon_list])

        output_file = f"{target.capitalize()}.GOA.Extracted.gaf"
        output_path = os.path.join(gaf_dir, output_file)

        try:
            with open(output_path, 'w') as outfile:
                # Uses zgrep as its faster than gzip parsing in python
                # pv utility tracks extraction process by tracking the bytes parsed in the input file
                cmd = f"pv -p -t -e -r {input_path} | zgrep -a '{taxon_pattern}'"
                subprocess.run(cmd, shell=True, stdout=outfile, check=True)

            print(f"{target.capitalize()} successfully extracted into {output_file}!\n")

        except subprocess.CalledProcessError as e:
            print(f"Extraction failed for {target.capitalize()}: {e}")

# FUNCTION: Download files from Xenbase: Xenbase GPI, genepage to gene id mapping, and ortholog mapping files
# Contains Xenbase URLs
def download_xenbase(args):
    print("Downloading Xenbase files...")

    if args.xen_download:
        # Download GPI file & unzip
        gpi_url = 'https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz'
        gpi_filepath = os.path.join(xb_dir, "Xenbase.gpi.gz")
        download_w_progress(gpi_url, gpi_filepath)
        unzip(gpi_filepath, remove_og=True)

    if args.ortho_download:
        # Download genepage to gene ID map
        # NOTE: line 11651 may be uncorrectly tabbed ("vma22  .L" instead of "vma22.L"); Check before using if redownloading
        genepage_to_gene_url = 'https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping_chd.txt'
        genepage_to_gene_filepath = os.path.join(xb_dir, "Xenbase_Genepage_To_GeneId.txt")
        download_w_progress(genepage_to_gene_url, genepage_to_gene_filepath)

        # Download genepage to ortholog NCBI ID map
        ortholog_map_url = 'https://xenbase-bio1.ucalgary.ca/cgi-bin/reports/genepage_entrez_orthologs.cgi'
        ortholog_map_filepath = os.path.join(ncbi_map_dir, "Xenopus_NCBI_Orthologs.tsv")
        download_w_progress(ortholog_map_url, ortholog_map_filepath)
        
        # Replace header with modified species names
        tmp_path = ortholog_map_filepath + ".tmp"
        with open(ortholog_map_filepath, 'r', encoding='utf-8') as infile, open(tmp_path, 'w', encoding='utf-8') as outfile:
            outfile.write("umbrella_id\ttxtrop\thuman\tmouse\trat\tzebrafish\tchicken\tdrosophila\tworm\n")
            lines = infile.readlines()      # Skip original header
            outfile.writelines(lines[1:])
        os.replace(tmp_path, ortholog_map_filepath)

    '''# Download current xenbase gaf & unzip:
    gaf_url = 'https://download.xenbase.org/xenbase/GenePageReports/xenbase.gaf.gz'
    gaf_filepath = os.path.join(xb_dir, "Xenbase.gaf.gz")
    download_w_progress(gaf_url, gaf_filepath)
    unzip(gaf_filepath, remove_og=True)'''

# FUNCTION: Download NCBI to uniprot mapping files for orthologs
# Contains uniprot URL
def download_maps():
    print("Downloading mapping files...")
    targets = [key for key in species_taxon_map.keys() if key != "xenopus"]

    for target in targets:
        # Map names to match partial species name used in uniprot urls (where required)
        species_map = {
            'chicken': 'chick',
            'zebrafish': 'danre',
            'drosophila': 'drome'
        }

        species_name = species_map.get(target, target).upper() # Default to original target value if target not present in species map
        species_taxon = species_taxon_map[target][0]  # Assuming one taxon per ortholog species
        
        #Build URL
        url = f'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{species_name}_{species_taxon}_idmapping.dat.gz'
        
        map_file = f"{target.capitalize()}_NCBI_Mapping.tsv"
        map_zipped = os.path.join(ncbi_map_dir, f"{map_file}.gz")
        # chunk size = <BYTE SIZE> (To modify number of bytes streamed in a chunk during download; default = 10*1024*1024 (10MB))

        download_w_progress(url, map_zipped)
        unzip(map_zipped, remove_og=True, line_filter=lambda line: "GeneID" in line)

# FUNCTION: Show progress while downloading
def download_w_progress(url, output_path, chunk_size=10*1024*1024):
    response = requests.get(url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get('content-length', 0))
    
    with open(output_path, 'wb') as f, tqdm(
        total=total_size,
        unit='B',
        unit_scale=True,
        unit_divisor=1024,
        desc=f"Downloading {os.path.basename(url)}",
        leave=True,
    ) as progress:
        for chunk in response.iter_content(chunk_size=chunk_size):
            f.write(chunk)
            progress.update(len(chunk))

    tqdm.write(f"Saved to {os.path.basename(output_path)}!\n")

# FUNCTION: Unzip .gz file; Removes .gz from filename if new filename not specified
def unzip(filepath, out_file=None, remove_og=False, line_filter=None):
    if not out_file:
        outpath = filepath.removesuffix(".gz")
    else:
        outpath = out_file

    mode = 'rt' if line_filter else 'rb'

    with gzip.open(filepath, mode) as in_file, open(outpath, 'w' if line_filter else 'wb') as out_file:
        if line_filter:
            for line in in_file:
                if line_filter(line):
                    out_file.write(line)
        else:
            shutil.copyfileobj(in_file, out_file)

    if remove_og:
        os.remove(filepath)
    print(f"Extracted to {outpath}!\n")
    return outpath

# FUNCTION: Adds annotations between 2 (GAF) files together
def combine_annotations(file_1, file_2, output_file, encoding="utf-8"):
    filter_func = lambda line: not line.startswith("!")
    temp1 = unzip(file_1, line_filter=filter_func)
    temp2 = unzip(file_2, line_filter=filter_func)

    try:
        with open(temp1, 'r', encoding=encoding) as f1, \
             open(temp2, 'r', encoding=encoding) as f2, \
             open(output_file, 'w', encoding=encoding) as out:
            out.writelines(f1.readlines())
            out.writelines(f2.readlines())
    finally:
        os.remove(temp1)
        os.remove(temp2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--xen_download', action='store_true')
    parser.add_argument('--ortho_download', action='store_true')
    parser.add_argument('--extract', nargs='?', const='all', help='Species name or "all"')
    parser.add_argument('--log', action='store_true')
    args = parser.parse_args()

    # Define folder paths (NOTE: If these are modified, must also change paths in goa_parsing.py)
    HOME = os.path.expanduser("~/xenbase-gaf-pipeline")
    input_dir = os.path.join(HOME, "input-files")
    gaf_dir = os.path.join(input_dir, "goa-gafs")
    ncbi_map_dir = os.path.join(input_dir, "ncbi-maps")
    xb_dir = os.path.join(input_dir, "xenbase-files")
    out_dir = os.path.join(HOME, "output-files")
    set_folders()

    # Redirect output to log file
    if args.log:
        log_path = os.path.join(out_dir, "script-logs/get_files.log")
        if args.xen_download:
            log = open(log_path, "wt")
        else:
            log = open(log_path, "at")
        sys.stdout = log
        sys.stderr = log    #tqdm & requests use sterr for progress bars

    # Set date and downloaded GOA filename
    print(f"Date & time of script execution: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")
    GOA_FILENAME = f"goa_uniprot_all.gaf.gz"

    # Set species to taxon map:
    species_taxon_map = {
        "xenopus": ["8355", "8364"],
        "human": ["9606"],
        "mouse": ["10090"],
        "rat": ["10116"],
        "chicken": ["9031"],
        "zebrafish": ["7955"],
        "drosophila": ["7227"]
    }

    # Download/extract files specified
    if args.xen_download or args.ortho_download:
        download_goa(args)
        download_xenbase(args)
        if args.ortho_download:
            extract_from_goa("all")
            download_maps(args)

    elif args.extract:
        extract_from_goa(args.extract)