import pandas as pd
import re
import os
import gzip
from datetime import datetime

HOME = os.path.expanduser("~")
target_dir = os.path.join(HOME, "xenbase-gaf-pipeline/input-files")
os.chdir(target_dir)

# Combines code from jupyter notebook into 4 general functions
def main():
  step_1()
  step_2()
  #step_3()
  #step_4()

# Filters human gaf to only include experimentally based evidence codes with xenopus orthologs
# Creates goa_human_with_orthologs_xenbase.tsv (used in step_2) (human gaf with mapped xenbase ids)
def step_1():
  # Define the allowed evidence codes
  # This is all experimentally based codes except for IPI
  allowed_evidence_codes = ["EXP", "IDA", "IMP", "IGI", "IEP", "ISS", "TAS", "IC"]

  def has_uniprotkb_accession(text):
    """Checks if the given text contains a UniProtKB accession pattern."""
    if pd.isnull(text) or text == "":  # Check for empty or missing values
      return False
    return bool(re.search(r"UniProtKB", text))

  # Load the human GO annotation file and filter by evidence codes
  goa_human = pd.read_csv("gafs/Human.GOA.Extracted.2025-07-18.gaf", sep="\t", header=None, comment="!", low_memory=False)
  goa_human = goa_human[goa_human[6].isin(allowed_evidence_codes)]  # Filter by evidence codes

  # Add the new filtering step here
  mask = goa_human[7].apply(has_uniprotkb_accession)  # Create boolean mask
  goa_human = goa_human[~mask]  # Keep rows where mask is False (no UniProtKB accession)

  # Extract the UniProtKB accessions from the GO annotation file
  human_uniprots = goa_human[1].unique()
  print(f"Number of unique UniProtKB accessions (Human): {len(human_uniprots)}")

  # Load the UniProt to GeneID mapping file
  idmapping = pd.read_csv("ncbi-map/Human_NCBI_Mapping.tsv", sep="\t", header=None)

  # Create a dictionary mapping UniProtKB accessions to GeneIDs
  uniprot_to_geneid = dict(zip(idmapping[0], idmapping[2]))
  print(f"Number of UniProtKB to GeneID mappings (Human): {len(uniprot_to_geneid)}")

  # Map UniProtKB accessions to GeneIDs in the human GO annotation file
  goa_human["GeneID"] = goa_human[1].map(uniprot_to_geneid)
  print(f"Number of GO annotations with GeneIDs (Human): {goa_human['GeneID'].notna().sum()}")

  # Load the Xenopus ortholog mapping file
  xenopus_orthologs = pd.read_csv("ncbi-map/Xenopus_NCBI_Orthologs.tsv", sep="\t")

  # Extract the human GeneIDs that have a Xenopus ortholog
  human_geneids_with_orthologs = xenopus_orthologs["human"].dropna().astype(int).unique()
  print(f"Number of human GeneIDs with Xenopus orthologs: {len(human_geneids_with_orthologs)}")

  # Filter the human GO annotation file to include only annotations with a Xenopus ortholog
  goa_human_with_orthologs = goa_human[goa_human["GeneID"].isin(human_geneids_with_orthologs)]
  #print(f"GOA Human with Orthologs shape: {goa_human_with_orthologs.shape}")

  # Save the filtered GO annotations to a file
  header = [
      "DB", "DB Object ID", "DB Object Symbol", "Qualifier", "GO ID",
      "DB:Reference", "Evidence Code", "With/From", "Aspect", "DB Object Name",
      "DB Object Synonym", "DB Object Type", "Taxon", "Date", "Assigned By",
      "Annotation Extension", "Gene Product ID"
  ]

  with open("goa_human_with_orthologs.tsv", "w") as f:
      f.write("\t".join(header) + "\n")
      goa_human_with_orthologs.to_csv(f, sep="\t", index=False, header=False)


  print("Human GO annotations with Xenopus orthologs saved to goa_human_with_orthologs.tsv")

  # Create a dictionary mapping human GeneIDs to Xenbase umbrella_ids
  geneid_to_xenbase = dict(zip(xenopus_orthologs["human"].fillna(-1).astype(int), xenopus_orthologs["umbrella_id"]))

  # Add a new column for Xenbase umbrella_id to the filtered GO annotations
  goa_human_with_orthologs["Xenbase_umbrella_id"] = goa_human_with_orthologs["GeneID"].map(geneid_to_xenbase)

  # Remove trailing zeros from "GeneID" column
  goa_human_with_orthologs["GeneID"] = goa_human_with_orthologs["GeneID"].astype(str).str.replace(r"\.0$", "", regex=True)

  # Save the filtered GO annotations with Xenbase IDs to a file, this has
  # additonal columns 18 and 19 containing the human Entrez ID and Xenbase
  # Genepage ID respectively.
  goa_human_with_orthologs.to_csv("goa_human_with_orthologs_xenbase.tsv", sep="\t", index=False, header=False)

  print("Human GO annotations with Xenopus orthologs and Xenbase IDs saved to goa_human_with_orthologs_xenbase.tsv")

# Builds x.trop gaf with human ortholog info
def step_2():
  # Constants
  PLACEHOLDER_REF = "GO_REF:PLACEHOLDER"
  TAXON_ID = "taxon:8364"
  EVIDENCE = "ISO"
  SOURCE = "Xenbase"
  DATE = datetime.now().strftime("%Y%m%d")

  # Load filtered annotations
  goa_human = pd.read_csv("goa_human_with_orthologs_xenbase.tsv", sep="\t", header=None)
  goa_human.columns = list(range(goa_human.shape[1]))

  print(f"✅ Loaded GOA with {goa_human.shape[0]} annotations")

  # Rename for readability
  goa_human["GeneID"] = goa_human[17].astype(str).str.replace(r"\.0$", "", regex=True)
  goa_human["Xenbase_umbrella"] = goa_human[18]

  goa_human["Xenbase_umbrella"] = "XB-GENEPAGE-" + goa_human[18].astype(str).str.strip()

  # Load genepage to Xenopus tropicalis gene mapping
  column_names = [
      "umbrella_id", "umbrella_symbol",
      "x_tropicalis_gene_id", "x_tropicalis_symbol",
      "x_laevis_L_gene_id", "x_laevis_L_symbol",
      "x_laevis_S_gene_id", "x_laevis_S_symbol"
  ]

  genepage_df = pd.read_csv(
      "Xenbase_Genepage_To_GeneId.txt",
      sep="\t",
      header=None,
      names=column_names,
      skiprows=1,
      dtype=str,
      engine="python"
  )

  print(f"✅ Loaded {genepage_df.shape[0]} Xenbase genepage mappings")

  # Build umbrella → X. tropicalis gene ID map
  genepage_to_xbgene = {
      row["umbrella_id"]: row["x_tropicalis_gene_id"]
      for _, row in genepage_df.iterrows()
      if pd.notnull(row["x_tropicalis_gene_id"]) and str(row["x_tropicalis_gene_id"]).startswith("XB-GENE")
  }

  print(f"✅ Mapped {len(genepage_to_xbgene)} umbrella → x.tropicalis IDs")

  # Load Xenbase GPI file (gene ID → symbol)
  xbgene_info = {}

  with open("Xenbase.gpi", "rt", encoding="latin-1") as f:
      count = 0
      for line in f:
          if line.startswith("!"):
              continue
          parts = line.strip().split("\t")
          if len(parts) < 3:
              continue
          db_object_id = parts[1]
          symbol = parts[2]
          xbgene_info[db_object_id] = {
              "symbol": symbol,
              "name": parts[3] if len(parts) > 3 else "",
              "synonym": parts[4] if len(parts) > 4 else "",
              "type": parts[5] if len(parts) > 5 else "protein"
          }
          count += 1

  print(f"✅ Loaded {count} gene records from GPI file")

  # Build new GAF lines
  gaf_lines = []
  skipped_missing_mappings = 0
  skipped_missing_gpi = 0
  mapped_successfully = 0

  for idx, row in goa_human.iterrows():
      umbrella_id = row["Xenbase_umbrella"]

      if pd.isna(umbrella_id) or umbrella_id not in genepage_to_xbgene:
          skipped_missing_mappings += 1
          continue

      xbgene = genepage_to_xbgene[umbrella_id]
      xbgene_curie = f"Xenbase:{xbgene}"

      gene_info = xbgene_info.get(xbgene)

      if row[2] == 'A1CF':
         print(f'Uniprot: {row[1]}, symbol: {row[2]}')
      if gene_info is None:
          skipped_missing_gpi += 1
          continue
      gaf_line = [
          "Xenbase",
          xbgene_curie,
          gene_info["symbol"],
          row[3],
          row[4],
          PLACEHOLDER_REF,
          EVIDENCE,
          f"UniProtKB:{row[1]}",
          row[8],
          gene_info["name"].strip(),
          gene_info["synonym"].strip(),
          gene_info["type"],
          TAXON_ID,
          DATE,
          SOURCE
      ]
      gaf_lines.append(tuple(gaf_line))
      mapped_successfully += 1

  print("Sample umbrella IDs from GOA file:")
  print(goa_human["Xenbase_umbrella"].dropna().unique()[:10])

  print(f"🔍 Mapped successfully: {mapped_successfully}")
  print(f"⚠️ Skipped due to missing x.tropicalis mapping: {skipped_missing_mappings}")
  print(f"⚠️ Skipped due to missing GPI data: {skipped_missing_gpi}")

  # Deduplicate
  df_gaf = pd.DataFrame(gaf_lines, columns=[
      "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference",
      "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "Synonym",
      "DB_Object_Type", "Taxon", "Date", "Assigned_By"
  ]).drop_duplicates(subset=["DB_Object_ID", "GO_ID", "With_From"])

  # Save to file
  df_gaf.to_csv("xenbase_from_human.gaf", sep="\t", header=False, index=False, quoting=3)
  print(f"✅ Xenopus GAF written with {len(df_gaf)} non-redundant entries.")

# Collapses gaf to genus level (replaces gene ids with umbrella ids)
def step_3():
  # Constants
  TAXON_GENEPAGE = "taxonid:262014"
  GAF_COLS = [
      "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference",
      "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "Synonym",
      "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
  ]

  # --- Load genepage-to-gene mapping ---
  column_names = [
      "umbrella_id", "umbrella_symbol",
      "x_tropicalis_gene_id", "x_tropicalis_symbol",
      "x_laevis_L_gene_id", "x_laevis_L_symbol",
      "x_laevis_S_gene_id", "x_laevis_S_symbol"
  ]
  genepage_df = pd.read_csv(
      "Xenbase_Genepage_To_GeneId.txt",
      sep="\t",
      header=None,
      names=column_names,
      skiprows=1,
      dtype=str
  )

  umbrella_id_to_symbol = {
      row["umbrella_id"]: row["umbrella_symbol"]
      for _, row in genepage_df.iterrows()
      if pd.notna(row["umbrella_id"]) and pd.notna(row["umbrella_symbol"])
  }
  # Invert all gene ID → umbrella mappings
  gene_to_umbrella = {}
  for _, row in genepage_df.iterrows():
      umbrella = row["umbrella_id"]
      for gene_col in ["x_tropicalis_gene_id", "x_laevis_L_gene_id", "x_laevis_S_gene_id"]:
          gene = row[gene_col]
          if pd.notna(gene) and gene.startswith("XB-GENE"):
              gene_to_umbrella[gene] = umbrella

  print(f"✅ Built mapping for {len(gene_to_umbrella)} gene IDs → umbrella IDs")

  # --- Load species-specific Xenbase GAF ---
  gaf = pd.read_csv(
      "gafs/Xenbase.gaf",
      sep="\t",
      comment="!",
      header=None,
      names=GAF_COLS,
      dtype=str
  )

  print(f"✅ Loaded {len(gaf)} Xenbase annotations")

  # --- Map each gene ID to umbrella ID ---
  # DB_Object_ID is like 'XB-GENE-12345678', so it maps directly
  gaf["Mapped_Umbrella"] = gaf["DB_Object_ID"].map(gene_to_umbrella)

  # Filter successfully mapped entries
  mapped = gaf[~gaf["Mapped_Umbrella"].isna()].copy()
  unmapped = gaf[gaf["Mapped_Umbrella"].isna()]
  print(f"🔍 Mapped {len(mapped)} annotations to umbrella IDs")
  print(f"⚠️ Skipped {len(unmapped)} with no umbrella mapping")

  # --- Update GAF fields to reflect umbrella ID context ---
  mapped["DB_Object_ID"] = mapped["Mapped_Umbrella"]
  mapped["DB_Object_Symbol"] = mapped["Mapped_Umbrella"].map(umbrella_id_to_symbol)
  mapped["Taxon"] = TAXON_GENEPAGE

  # Optional: also update DB_Object_Symbol or Assigned_By if needed

  # --- Deduplicate based on umbrella ID and GO term ---
  df_umbrella_gaf = mapped[GAF_COLS].drop_duplicates(subset=["DB_Object_ID", "GO_ID"])

  # --- Save output ---
  df_umbrella_gaf.to_csv("xenbase_species_collapsed_to_genepage.gaf", sep="\t", header=False, index=False, quoting=3)
  print(f"✅ Collapsed species-specific GAF written with {len(df_umbrella_gaf)} umbrella-level annotations.")
  
  # --- Debug samples to check why mapping fails ---
  print("\n🔎 Sample DB_Object_IDs from GAF:")
  print(gaf["DB_Object_ID"].dropna().unique()[:10])

  print("\n🔎 Sample gene_to_umbrella keys:")
  print(list(gene_to_umbrella.keys())[:10])

  # Check overlap directly
  gaf_ids = set(gaf["DB_Object_ID"].dropna().unique())
  mapping_keys = set(gene_to_umbrella.keys())

  overlap = gaf_ids.intersection(mapping_keys)
  print(f"\n🔁 Overlapping gene IDs: {len(overlap)}")
  if overlap:
      print("✅ Sample overlap:", list(overlap)[:5])
  else:
      print("❌ No overlap found – mismatch likely.")

# Merge gafs produced in step 2 & 3 (updating to use umbrella id) for final output
def step_4():
  # Constants
  TAXON_GENEPAGE = "taxonid:262014"
  GAF_COLS = [
      "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference",
      "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "Synonym",
      "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
  ]

  # --- Load genepage-to-gene mapping ---
  column_names = [
      "umbrella_id", "umbrella_symbol",
      "x_tropicalis_gene_id", "x_tropicalis_symbol",
      "x_laevis_L_gene_id", "x_laevis_L_symbol",
      "x_laevis_S_gene_id", "x_laevis_S_symbol"
  ]
  genepage_df = pd.read_csv(
      "Xenbase_Genepage_To_GeneId.txt",
      sep="\t",
      header=None,
      names=column_names,
      skiprows=1,
      dtype=str
  )

  # Build gene → umbrella + symbol maps
  gene_to_umbrella = {}
  umbrella_symbol_map = {}

  for _, row in genepage_df.iterrows():
      umbrella = row["umbrella_id"]
      symbol = row["umbrella_symbol"]
      for gene_col in ["x_tropicalis_gene_id", "x_laevis_L_gene_id", "x_laevis_S_gene_id"]:
          gene = row[gene_col]
          if pd.notna(gene) and gene.startswith("XB-GENE"):
              gene_to_umbrella[gene] = umbrella
              umbrella_symbol_map[umbrella] = symbol

  print(f"✅ Loaded mapping for {len(gene_to_umbrella)} gene IDs to umbrella IDs")

  # --- Load both GAFs ---
  gaf_human = pd.read_csv("xenbase_from_human.gaf", sep="\t", header=None, names=GAF_COLS, dtype=str)
  gaf_species = pd.read_csv("xenbase_species_collapsed_to_genepage.gaf", sep="\t", header=None, names=GAF_COLS, dtype=str)

  print(f"✅ Loaded {len(gaf_human)} human-derived annotations")
  print(f"✅ Loaded {len(gaf_species)} species-specific umbrella annotations")

  # --- Normalize human GAF to umbrella IDs ---
  gaf_human["Original_Gene_ID"] = gaf_human["DB_Object_ID"].str.replace("Xenbase:", "", regex=False)
  # Map to umbrella
  gaf_human["Mapped_Umbrella"] = gaf_human["Original_Gene_ID"].map(gene_to_umbrella)
  gaf_human["Mapped_Symbol"] = gaf_human["Mapped_Umbrella"].map(umbrella_symbol_map)
  print(gaf_human.columns)

  # Filter only successfully mapped entries
  mapped_human = gaf_human[~gaf_human["Mapped_Umbrella"].isna()].copy()

  # Apply umbrella-level changes
  mapped_human["DB_Object_ID"] = mapped_human["Mapped_Umbrella"]
  mapped_human["DB_Object_Symbol"] = mapped_human["Mapped_Symbol"]
  mapped_human["Taxon"] = TAXON_GENEPAGE

  print(f"🔁 Successfully remapped {len(mapped_human)} human annotations to umbrella IDs")

  # --- Merge and deduplicate ---
  combined = pd.concat([mapped_human[GAF_COLS], gaf_species], ignore_index=True)
  combined = combined.drop_duplicates(subset=[
      "DB_Object_ID", "GO_ID", "Evidence_Code", "Qualifier", "With_From"
  ])
  # --- Write final GAF ---
  combined.to_csv("xenbase_final_genepage_merged.gaf", sep="\t", header=False, index=False, quoting=3)
  print(f"✅ Final merged genepage-level GAF written with {len(combined)} annotations.")

main()