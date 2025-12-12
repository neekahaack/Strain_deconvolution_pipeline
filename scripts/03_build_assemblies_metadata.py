import pandas as pd
from pathlib import Path
import sys

ASM_DIR = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assembled_genomes")
CSV_PATH = Path("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/41587_2018_9_MOESM3_ESM.csv")
META_OUT = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assemblies_metadata.tsv")

print("Checking inputs...")

#Check if assembled genomes are present
if not ASM_DIR.exists():
    print(f"ERROR: Assembly directory does not exist: {ASM_DIR}")
    sys.exit(1)

if not CSV_PATH.exists():
    print(f"ERROR: Metadata CSV not found: {CSV_PATH}")
    sys.exit(1)

#Load in CSV to df variable
print("Loading CSV...")
df = pd.read_csv(CSV_PATH, sep=";", skiprows=1)

if "Identifier" not in df.columns:
    print("ERROR: Could not find 'Identifier' column in CSV.")
    print("Columns:", df.columns)
    sys.exit(1)

#Rename Identifier to sample_ID to match 
df = df.rename(columns={"Identifier": "sample_id"})

#Filter HBC only
df = df[df.iloc[:, 1] == "HBC"]

#Initialized records and missing variables
records = []
missing = []

print("Processing assemblies...")

for sample_dir in ASM_DIR.iterdir():
    contigs = sample_dir / "contigs.fasta"

    #Check assemblies to make sure theyre useable
    if not contigs.exists() or contigs.stat().st_size == 0:
        missing.append(sample_dir.name)
        continue

    sample_id = sample_dir.name     #the sample id is the directory name
    row = df[df["sample_id"] == sample_id]  #filters for that to be true ^

    if len(row) == 0:
        species = "UNKNOWN"
    else:
        species = " ".join(str(row.iloc[0]["Sample"]).split()[:-1])

    records.append({
        "Sample_id": sample_id,
        "Species": species,
        "Assembly_path": str(contigs)
    })

meta_df = pd.DataFrame(records)
meta_df.to_csv(META_OUT, sep="\t", index=False)

print(f"\nMetadata written to: {META_OUT}")
print(f"Total assemblies recorded: {len(records)}")

if missing:
    print("\nMissing or broken assemblies:")
    for m in missing:
        print(" -", m)
else:
    print("\nAll assemblies present!")