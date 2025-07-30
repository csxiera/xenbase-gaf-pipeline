import argparse
import requests
import os
import subprocess
import gzip
import shutil
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

# NOTE: URL strings may change in the future! Update EBI/uniprot/xenbase URL string in appropriate function if/when this happens

# Expected download times (from July 2025):
#   - goa_uniprot_all.gaf.gz (~4hrs, 21GB)
#   - species-specific gafs (~18 min each)
#   - remaining files (>1 min each)

# FUNCTION: Make input folders if they dont yet exist
def set_folders():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(GAF_DIR, exist_ok=True)
    os.makedirs(NCBI_MAP_DIR, exist_ok=True)
    os.makedirs(XB_DIR, exist_ok=True)

# FUNCTION: Download GOA file with data for all species
# Contains EBI URL
def download_goa():
    print("Downloading GOA file from EBI...")
    url = f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
    output_path = os.path.join(GAF_DIR, "original-goa", GOA_FILENAME)
    # chunk size = <BYTE SIZE> (To modify number of bytes streamed in a chunk during download; default = 10*1024*1024 (10MB))

    download_w_progress(url, output_path=output_path)

# FUNCTION: Extract GAF data for each species into seperate files
def extract_from_goa(species):
    input_path = os.path.join(GAF_DIR, "original-goa", GOA_FILENAME)
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

        output_file = f"{target.capitalize()}.GOA.Extracted.{DATE}.gaf"
        output_path = os.path.join(GAF_DIR, output_file)

        try:
            with open(output_path, 'w') as outfile:
                '''subprocess.run(["zgrep", "-a", taxon_pattern, input_path],
                            stdout=outfile,
                            check=True)'''
                
                # Uses zgrep as its faster than gzip parsing in python
                # pv utility tracks extraction process by tracking the bytes parsed in the input file
                # FIX!!: Remove path to pv once installed by kamran
                cmd = f"export PATH=$HOME/.local/bin:$PATH && pv -p -t -e -r {input_path} | zgrep -a '{taxon_pattern}'"
                subprocess.run(cmd, shell=True, stdout=outfile, check=True)

            print(f"{target.capitalize()} successfully extracted into {output_file}!\n")

        except subprocess.CalledProcessError as e:
            print(f"Extraction failed for {target.capitalize()}: {e}")

# FUNCTION: Download NCBI to uniprot mapping files
# Contains uniprot URL
def download_maps():
    print("Downloading mapping files...\n")
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
        
        output_file = f"{target.capitalize()}_NCBI_Mapping.tsv"
        output_zipped = os.path.join(NCBI_MAP_DIR, f"{output_file}.gz")
        # chunk size = <BYTE SIZE> (To modify number of bytes streamed in a chunk during download; default = 10*1024*1024 (10MB))

        download_w_progress(url, output_zipped)

        with gzip.open(output_zipped, 'rt') as gz_file, open(output_zipped.removesuffix(".gz"), 'w') as out_file:
            for line in gz_file:
                if 'GeneID' in line:
                    out_file.write(line)
        os.remove(output_zipped)
        print(f"Extracted to {output_file}!\n")


# FUNCTION: Download Xenbase GPI, genepage to gene id mapping, and ortholog mapping files
# Contains Xenbase URLs
def download_xenbase():
    print("Downloading Xenbase files...\n")

    # Download GPI file
    gpi_url = 'https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz'
    gpi_output = os.path.join(XB_DIR, "Xenbase.gpi.gz")
    download_w_progress(gpi_url, gpi_output)
    
    # Unzip GPI and remove .gz version
    with gzip.open(gpi_output, 'rb') as gz_file, open(gpi_output.removesuffix(".gz"), 'wb') as out_file:
        shutil.copyfileobj(gz_file, out_file)
    os.remove(gpi_output)
    print(f"Extracted to {gpi_output.removesuffix(".gz")}!\n")

    # Download genepage to gene ID map
    # NOTE: line 11651 may be uncorrectly tabbed ("vma22  .L" instead of "vma22.L"); Check before using if redownloading
    genepage_to_gene_url = 'https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping_chd.txt'
    genepage_to_gene_output = os.path.join(XB_DIR, "Xenbase_Genepage_To_GeneId.txt")
    download_w_progress(genepage_to_gene_url, genepage_to_gene_output)

    # Download genepage to ortholog NCBI ID map
    ortholog_map_url = 'https://xenbase-bio1.ucalgary.ca/cgi-bin/reports/genepage_entrez_orthologs.cgi'
    ortholog_map_output = os.path.join(NCBI_MAP_DIR, "Xenopus_NCBI_Orthologs-test.tsv")
    download_w_progress(ortholog_map_url, ortholog_map_output)
    
    # Replace header with modified species names
    tmp_path = ortholog_map_output + ".tmp"
    with open(ortholog_map_output, 'r', encoding='utf-8') as infile, open(tmp_path, 'w', encoding='utf-8') as outfile:
        outfile.write("umbrella_id\ttxtrop\thuman\tmouse\trat\tzebrafish\tchicken\tdrosophila\tworm\n")
        lines = infile.readlines()      # Skip original header
        outfile.writelines(lines[1:])
    os.replace(tmp_path, ortholog_map_output)

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--goa', action='store_true')
    parser.add_argument('--extract', nargs='?', const='all', help='Species name or "all"')
    parser.add_argument('--maps', action='store_true')
    parser.add_argument('--xenbase', action='store_true')
    parser.add_argument('--date', help='Optional date in YYYY-MM-DD format (Use if regenerating species-specific GAFs from older GOA download)')
    args = parser.parse_args()

    # Define file paths (NOTE: If these are modified, must also change paths in goa_parsing.py)
    HOME = Path.home()
    DATA_DIR = f"{HOME}/xenbase-gaf-pipeline/input-files"
    GAF_DIR = f"{DATA_DIR}/goa-gafs"
    NCBI_MAP_DIR = f"{DATA_DIR}/ncbi-maps"
    XB_DIR = f"{DATA_DIR}/xenbase-files"
    set_folders()

    # Set date and downloaded GOA filename
    DATE = args.date if args.date else datetime.today().strftime('%Y-%m-%d')
    GOA_FILENAME = f"goa_uniprot_all.{DATE}.gaf.gz"

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
    if args.goa: download_goa()
    if args.extract: extract_from_goa(args.extract)
    if args.xenbase: download_xenbase()
    if args.maps: download_maps()
