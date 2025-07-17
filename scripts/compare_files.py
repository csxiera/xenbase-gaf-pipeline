def read_lines(filepath):
    with open(filepath, 'r', encoding='latin1') as f:
        return set(line.rstrip('\n') for line in f)

def compare_files(file_a, file_b):
    output_a=f'only_{file_a}.txt'
    output_b=f'only_{file_b}.txt'

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


goa_file = 'Xenopus.GOA.2025-06-24.gaf'
matched_file = f'matched_annotations.{goa_file}'
ebi_only_file = 'xenbase.EBI.only.2.2.gaf'

compare_files(matched_file, ebi_only_file)