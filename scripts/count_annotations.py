import os, csv

HOME = os.path.expanduser("~/xenbase-gaf-pipeline")
input_dir = os.path.join(HOME, "output-files/ortho-gafs")

#file_name = "Xenopus.GOA.Extracted.2025-07-18.gaf"
#file_name = "Human.GOA.Extracted.2025-07-18.gaf"
#file_name = "Master_Orthologs.gaf"
file_name = "Xenbase_from_Human.gaf"

input_file = os.path.join(input_dir, file_name)
output_file = os.path.join(HOME, "output-files/count_annotations.tmp")

annotation_count = 0
count_on_column = "all"
#count_on_column = "evidence"
#count_by_values = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]
count_by_values = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "TAS", "IC"]

with open(input_file, 'r', encoding='latin1') as f_in, \
    open(output_file, 'w', encoding='latin1', newline='') as f_out:
    
    reader = csv.reader(f_in, delimiter='\t')
    writer = csv.writer(f_out, delimiter='\t')
    writer.writerow(f"! Annotation file to count: {input_file}\n")
    if count_on_column != "all"
        writer.writerow(f"! Counting annotations in {count_on_column} column with values matching one of the following: [{','.join(count_by_values)}]")
    else:
        writer.writerow(f"! Counting all annotations")

    for fields in reader:
        if not fields or fields[0].startswith('!') or len(fields) < 17:
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

        if count_on_column != "all"
            if records[count_on_column] in count_by_values:
                annotation_count += 1
                writer.writerow(fields)
        
        else:
            annotation_count += 1
        
        
print(f'There are {annotation_count} annotations with experimentally supported evidence in {input_file}')
