import subprocess
from pathlib import Path
import pandas as pd
import requests

#Where to store FASTQ files
fastq_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/raw_reads")
fastq_dir.mkdir(parents=True, exist_ok=True)

#Load the original CSV (1354 rows)
df = pd.read_csv(
    "/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/41587_2018_9_MOESM3_ESM.csv",
    sep=";",
    skiprows=1
)

#We assume the column containing run accessions is named "Identifier"
accessions = df["Identifier"].dropna().astype(str)

print(f"Found {len(accessions)} accessions in the CSV.")

for run in accessions:
    run = run.strip()
    if not run:
        continue

    print(f"\n=== Getting FASTQ links for {run} ===")

    #Query ENA for FASTQ URLs
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={run}&result=read_run&fields=fastq_ftp&format=tsv&download=false"
    )

    try:
        r = requests.get(url)
    except Exception as e:
        print(f"Error querying ENA: {e}")
        continue

    lines = r.text.strip().split("\n")
    if len(lines) < 2 or not lines[1].strip():
        print("No FASTQ found for this accession.")
        continue

    ftp_field = lines[1].split("\t")[-1]
    ftp_links = [f"ftp://{x}" for x in ftp_field.split(";") if x.strip()]

    if not ftp_links:
        print("No FASTQ links in ENA response.")
        continue

    #Download each FASTQ file
    for link in ftp_links:
        filename = fastq_dir / link.split("/")[-1]

        if filename.exists() and filename.stat().st_size > 0:
            print(f"  Already have {filename.name}, skipping.")
            continue

        print(f"  Downloading {filename.name} ...")
        subprocess.run(["wget", "-q", "-O", str(filename), link])

print("\n Done downloading FASTQ files for all species!")

# ================================================================
# Check for missing isolates based on the CSV
# ================================================================

print("\n Checking for missing isolates...")

#expected list from CSV (HBC only)
hbc_df = df[df.iloc[:,1] == "HBC"]      # column 2 is index 1
expected_ids = set(hbc_df["Identifier"].dropna().astype(str).str.strip())

#what was actually downloaded (_1.fastq.gz)
downloaded_ids = {
    p.name.replace("_1.fastq.gz", "")
    for p in fastq_dir.glob("*_1.fastq.gz")
}

missing = sorted(expected_ids - downloaded_ids)

if missing:
    print(f"⚠ Missing isolates: {len(missing)}")
    for m in missing:
        print("  -", m)
else:
    print("No isolates missing!")

print("\nScript complete.")