import subprocess
from pathlib import Path
import pandas as pd

# ---------------------------
# PATHS
# ---------------------------
input_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/raw_reads")
output_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assembled_genomes")
output_dir.mkdir(exist_ok=True)

# Path to the original metadata CSV (1354 isolates)
csv_path = Path("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/41587_2018_9_MOESM3_ESM.csv")

# Track successful assemblies
assembled_samples = []

# ---------------------------
# ASSEMBLE GENOMES
# ---------------------------
for r1 in sorted(input_dir.glob("*_1.fastq.gz")):
    sample = r1.name.replace("_1.fastq.gz", "")
    r2 = input_dir / f"{sample}_2.fastq.gz"

    if not r2.exists():
        print(f"Missing reverse read for {sample}, skipping.")
        continue

    sample_out = output_dir / sample
    contigs = sample_out / "contigs.fasta"

    # If contigs exist, this isolate is already assembled
    if contigs.exists():
        print(f"{sample} already assembled, skipping.")
        assembled_samples.append(sample)
        continue

    print(f"Assembling {sample}...")

    cmd = [
        "spades.py",
        "-1", str(r1),
        "-2", str(r2),
        "-o", str(sample_out),
        "--threads", "8",
        "--memory", "64",
        "--isolate"
    ]

    result = subprocess.run(cmd)

    # Only count as assembled if contigs.fasta exists
    if contigs.exists():
        print(f"Finished assembling {sample}.")
        assembled_samples.append(sample)
    else:
        print(f"❌ Assembly failed for {sample} (no contigs.fasta created).")

print("\nAll assemblies complete!")

# ---------------------------
# BUILD METADATA TABLE
# ---------------------------
print("\nBuilding assemblies_metadata.tsv ...")

df = pd.read_csv(csv_path, sep=";", skiprows=1)
df = df.rename(columns={"Identifier": "sample_id"})

# ⭐ FILTER HBC ONLY (important)
df = df[df.iloc[:, 1] == "HBC"]

# Detect species column
species_col = None
for c in df.columns:
    if "species" in c.lower() or "organism" in c.lower():
        species_col = c
        break

if species_col is None:
    raise ValueError("Could not detect species column in CSV. Inspect df.columns.")

records = []

for sample_id in assembled_samples:
    row = df[df["sample_id"] == sample_id]

    if len(row) == 0:
        species = "UNKNOWN"
    else:
        species = row[species_col].iloc[0]

    assembly_path = output_dir / sample_id / "contigs.fasta"

    if assembly_path.exists():
        records.append({
            "sample_id": sample_id,
            "species": species,
            "assembly_path": str(assembly_path)
        })

meta_df = pd.DataFrame(records)
meta_out = output_dir.parent / "assemblies_metadata.tsv"
meta_df.to_csv(meta_out, sep="\t", index=False)

print(f"Metadata saved to: {meta_out}")
print(f"Total assembled genomes recorded: {len(meta_df)}")

