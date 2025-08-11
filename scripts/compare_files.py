import os
from pathlib import Path

def read_lines(filepath, columns=None):
    lines = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!'):
                continue
            full_line = line.rstrip('\n')
            if columns is not None:
                fields = full_line.split('\t')
                key = tuple(fields[i] for i in columns if i < len(fields))
            else:
                key = full_line
            lines[key] = full_line
    return lines

def compare_files(file_a, file_b, out_dir, columns=None):
    print(f'\nComparing lines between files:')
    print(f'File A: {file_a}')
    print(f'File B: {file_b}\n')

    filename_a, _ = os.path.splitext(os.path.basename(file_a))      # Get file basename and strip extension
    filename_b, _ = os.path.splitext(os.path.basename(file_b))
    out_a = os.path.join(out_dir, f'{filename_a}_A_Only.tmp')
    out_b = os.path.join(out_dir, f'{filename_b}_B_Only.tmp')

    lines_a = read_lines(file_a, columns=columns)
    lines_b = read_lines(file_b, columns=columns)

    keys_a = set(lines_a.keys())
    keys_b = set(lines_b.keys())

    only_a_keys = keys_a - keys_b
    only_b_keys = keys_b - keys_a

    if only_a_keys:
        print(f'{len(only_a_keys)} unique entries found in {file_a}\n')
        with open(out_a, 'w', encoding='utf-8') as f:
            f.write(f'! Entries only in {file_a} (File A)\n')
            f.write(f'! -> Compared to {file_b} (File B)\n\n')
            for key in sorted(only_a_keys):
                f.write(lines_a[key] + '\n')
    else:
        print(f'There are no unique entries in {file_a}')

    if only_b_keys:
        print(f'{len(only_b_keys)} unique entries found in {file_b}\n')
        with open(out_b, 'w', encoding='utf-8') as f:
            f.write(f'! Entries only in {file_b} (File B)\n')
            f.write(f'! -> Compared to {file_a} (File A)\n\n')
            for key in sorted(only_b_keys):
                f.write(lines_b[key] + '\n')
    else:
        print(f'There are no unique entries in {file_b}')

def run_compare(preset=None, file_a=None, file_b=None, columns=None):
    out_dir = os.path.join(HOME, 'xenbase-gaf-pipeline/output-files/analysis')

    if preset == "xenbase-vs-new":
        dir1 = os.path.join(HOME, 'xenbase-gaf-pipeline/input-files/xenbase-files')
        dir2 = os.path.join(HOME, 'xenbase-gaf-pipeline/output-files')
        file_a = os.path.join(dir1, "Xenbase.gaf")
        file_b = os.path.join(dir2, "Xenbase.gaf")

    if not (file_a and file_b):
        raise ValueError("You must provide either a preset or both input filepaths")

    compare_files(file_a, file_b, out_dir, columns)

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
    #
    # SCRIPT FUNCTION: Compare 2 files (Usually GAF), write "unique" lines in each to seperate files named {file_name}_A/B_Only.tmp
    #
    # 1. Basic file comparison:
    #    - By default, annotations in each file will be compared on the whole line
    #
    # 2. File comparison on specific columns
    #    - Annotations in each file are compared based on their values in a given set of columns
    #    - Ex. If an annotation from each file matches on their values in columns 2,3, and 6, 
    #       regardless of whether they differ in other column values, they will be treated as matches
    #
    # Presets (optional shortcut flags for common files used):
    #
    #   preset="xenbase-vs-new"     -> Compare 'original' Xenbase.gaf (version currently available on site) to new gaf produced by pipeline
    #   preset="analysis"           -> Compare 2 files in analysis folder; file names required
    #   NOTE: If no preset is used, you must provide both file names and their directories
    #
    # NOTE: All output files (if created) will be saved to output-files/analysis
    
    HOME = Path.home()
    ANALYSIS_DIR = os.path.join(HOME, "xenbase-gaf-pipeline/output-files/analysis")
    os.makedirs(ANALYSIS_DIR, exist_ok=True)

    # Use to compare 'old' and 'new' xenbase gafs
    #run_compare(preset="xenbase-vs-new", columns=[1,4,5])       # Object ID, GO ID, Evidence Code

    # Use to compare non-preset files:
    filepath_a = os.path.join(HOME, "xenbase-gaf-pipeline/output-files/xenbase.EBI.only.2.2.gaf")
    filepath_b = os.path.join(HOME, "xenbase-gaf-pipeline/output-files/original-pipeline-output/xenbase.EBI.only.2.2.gaf")

    run_compare(file_a=filepath_a, file_b=filepath_b, columns=[0,1,2,3,4,5,6])

    # Use to clean up {...}_Only_A/B.tmp files
    #remove_tmp_files(output_dir)
