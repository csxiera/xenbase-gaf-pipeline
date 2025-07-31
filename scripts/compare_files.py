import os
from pathlib import Path

def read_lines(filepath):
    with open(filepath, 'r', encoding='latin1') as f:
        return set(line.rstrip('\n') for line in f)

def compare_files(file_a, file_b):
    output_a=f'{file_a}_only.txt'
    output_b=f'{file_b}_only.txt'

    lines_a = read_lines(file_a)
    lines_b = read_lines(file_b)

    only_a = lines_a - lines_b
    only_b = lines_b - lines_a

    with open(output_a, 'w', encoding='utf-8') as f:
        for line in sorted(only_a):
            f.write(line + '\n')

    with open(output_b, 'w', encoding='utf-8') as f:
        for line in sorted(only_b):
            f.write(line + '\n')

HOME = Path.home()
DIR = f"{HOME}/xenbase-gaf-pipeline/input-files/goa-gafs"

file1 = os.path.join(DIR, 'Human.GOA.Extracted.2025-07-18.gaf')
file2 = os.path.join(DIR, 'Human.GOA.Extracted.2025-07-30.gaf')

compare_files(file1, file2)