import os, csv
from pathlib import Path

def count(input_file, filter_column, filter_values=None, find_unique_values=False, find_in_column=None, output_file=None):
    annotation_count = 0
    unique_values = set()

    with open(input_file, 'r', encoding='latin1') as f_in:
        reader = csv.reader(f_in, delimiter='\t')
        writer = None

        print(f"Annotation file to count: {input_file}\n")
        # Only write to output file if counting by specific column/values (ie. create filtered file)
        if output_file:
            if not filter_column or not filter_values:
                print("No filtering cloumn/values specified. No output file will be produced")
                output_file = None
            else:
                print(f"Matching annotations will be written to {output_file}\n")
                f_out = open(output_file, 'w', encoding='latin1', newline='')
                writer = csv.writer(f_out, delimiter='\t')
                writer.writerow([f"! Annotation file to filter: {input_file}"])
                writer.writerow([f"! Filtering annotations in {filter_column} column with values matching one of the following: [{','.join(filter_values)}]"])
        elif find_unique_values:
            if not find_in_column and not filter_column:
                print('Error: Must provide column to find unique values in')
                return
        elif not filter_column:
                print("Counting all annotations, no output file will be produced\n")

        for fields in reader:
            if not fields or fields[0].startswith('!') or len(fields) < 17:
                continue

            target_col = find_in_column if find_in_column else filter_column

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

            if filter_column and filter_values:
                if records[filter_column] in filter_values:
                    annotation_count += 1
                    if find_unique_values and target_col in records:
                        unique_values.add(records[target_col])
                    if writer:
                        writer.writerow(fields)
            else:
                annotation_count += 1
                if find_unique_values and target_col in records:
                    unique_values.add(records[target_col])

        if output_file:
            f_out.close()

    if filter_column and filter_values:
        print(f'There are {annotation_count} annotations with "{filter_column}" values in {filter_values} for {input_file}')
    else:
        print(f'There are {annotation_count} annotations in {input_file}')

    if find_unique_values:
        print(f'\nThere are {len(unique_values)} unique values in the "{target_col}" column')
        if len(unique_values) < 50:
            print(f'Unique {target_col} values:')
            print(unique_values)

def run_count(preset=None, species=None, date=None, file_name=None, filter_column=None, filter_values=None, find_unique=False, find_in_column=None, output_file=None):
    output_dir = os.path.join(HOME, "xenbase-gaf-pipeline/output-files/analysis")

    if preset:
        species = species.capitalize() if species else None
        if preset == "input-goa":
            if not (species and date):
                raise ValueError("Both species and date must be provided for 'input-goa'")

            input_dir = os.path.join(HOME, "xenbase-gaf-pipeline/input-files/goa-gafs")
            file_name = f"{species}.GOA.Extracted.{date}.gaf"

            if not os.path.exists(input_file):
                raise FileNotFoundError(f"Input GAF file not found: {input_file}")

            output_file = os.path.join(output_dir, f"Count_Annotations.{species}.GOA.tmp")

        elif preset == "xenbase-gaf":
            input_dir = os.path.join(HOME, "xenbase-gaf-pipeline/input-files/xenbase-files")
            file_name = "Xenbase.gaf"
            output_file = os.path.join(output_dir, "Count_Annotations.Xenbase_Original_GAF.tmp")

        elif preset == "output-gaf":
            if not species:
                raise ValueError("Species must be provided for 'output-gaf'")

            input_dir = output_dir
            file_name = f"Xenbase_w_{species}.gaf" if species != "Xenopus" else "Xenbase.gaf"
            output_file = os.path.join(output_dir, f"Count_Annotations.{file_name.removesuffix(".gaf")}.tmp")

        elif preset == "analysis" and file_name:
            input_dir = output_dir
            if output_file:
                output_file = os.path.join(output_dir, output_file)

        else:
            raise ValueError(f"Unknown file preset type: '{preset}'")

        input_file = os.path.join(input_dir, file_name)

    # Set input and output files for non-preset option
    elif os.path.isabs(file_name):
        input_file = file_name
        if output_file:
            output_file = os.path.join(output_dir, os.path.basename(output_file)) # ensures correct output directory used even if full path is provided

    else:
        raise ValueError("You must provide either a preset and/or a file/filepath")
    
    count(input_file=input_file, filter_column=filter_column, filter_values=filter_values, find_unique_values=find_unique, find_in_column=find_in_column, output_file=output_file)

def remove_tmp_files(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith(".tmp"):
            file_path = os.path.join(folder_path, filename)
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

if __name__ == "__main__":
    # AUTHOR: C. Lenz
    # DATE LAST UPDATED: 2025-08-05
    #
    # SCRIPT FUNCTION: Filtering and counting of annotations in GAF files
    #
    # 1. Basic annotation counting:
    #    - By default, all annotations in the input file will be counted with no filtering applied
    #    - No ouput file created (would be identical to input)
    #
    # 2. Filtering annotations by column values:
    #    - Use 'filter_column' and 'filter_values' to count only rows where the specified column contains one of the given values
    #    - Optional: provide 'output_file' to write the filtered annotations to a new file
    #
    # 3. Finding unique values in a column:
    #    - Set 'find_unique=True' to find all unique values in a specified column
    #    - With filtering: set 'find_in_column' to specify which column to find unique values (within the filtered subset)
    #    - Without filtering: set either 'find_in_column' or 'filter_column" to specify which column to find unique values
    #    NOTE: If fewer than 50 unique values are found, they are printed to the terminal, otherwise only the count of unique values is printed
    #
    # Presets (optional shortcut flags for common files used):
    #
    #   preset="input-goa"     -> Use a species-specific GOA GAF (requires 'species' and 'date' arguments)
    #   preset="xenbase-gaf"   -> Use original Xenbase GAF (downloaded from site)
    #   preset="output-gaf"    -> Use a GAF file produced by Xenbase GAF pipeline (requires 'species' argument)
    #   preset="analysis"      -> Use a specific file in the analysis folder (requires 'file_name' argument)
    #   NOTE: If no preset is used, you must provide full file path as 'file_name'
    #
    # NOTE: All output files (if created) will be saved to output-files/analysis

    HOME = Path.home()
    output_dir = os.path.join(HOME, "xenbase-gaf-pipeline/output-files/analysis")
    os.makedirs(output_dir, exist_ok=True)

    # Use to analyze original Xenbase GAF:
    evidence_codes = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "TAS", "IC"]
    #run_count(preset="xenbase-gaf", filter_column="evidence", filter_values=evidence_codes)

    # Use for files in analysis folder (xenbase-gaf-pipeline/output-files/analysis):
    input_file = "Xenbase_A_Only.tmp"
    output_file = "XB-Genes_Only.tmp"
    evidence_codes=["IDA","IMP","ISS","ISO","TAS","NAS","IEA"]
    #run_count(preset="analysis", file_name=input_file, filter_column="evidence", filter_values="NAS", find_unique=True, find_in_column='assigned_by')

    input_file = "Xenbase.EBI.only.NEW.2025-07-18_B_Only.tmp"
    run_count(preset="analysis", file_name=input_file, find_unique=True, find_in_column="go_id")

    # Use for non-preset files (requires full file path):
    input_file = os.path.join(HOME, "GOA_pipeline/xenbase.EBI.only.2.2.gaf")
    #run_count(file_name=input_file)

    # Use to clean up Count_Annotations{...}.tmp files
    #remove_tmp_files(output_dir)
