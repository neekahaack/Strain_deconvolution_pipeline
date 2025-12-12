import os
import time
import urllib.request
from Bio import Entrez
import pandas as pd

# Set email, required by NCBI Entrez
Entrez.email = "neeka.haack@embl.de"

# Map project → input file
INPUT_FILES = {
    "PRJNA482748": "/g/typas/Personal_Folders/Neeka/Model/data/PRJNA482748_AssemblyDetails.txt",
    "PRJNA903559": "/g/typas/Personal_Folders/Neeka/Model/data/PRJNA903559_AssemblyDetails.txt",
}
BASE_OUTDIR = "/g/typas/Personal_Folders/Neeka/Model/data/assemblies/"  # parent folder

def fetch_fasta_from_assembly(assembly_id: str, out_dir: str):
    """Download genomic FASTA for one assembly accession into out_dir."""
    search = Entrez.esearch(db="assembly", term=assembly_id, retmode="xml")
    search_record = Entrez.read(search)
    uid_list = search_record.get("IdList", [])
    if not uid_list:
        print(f"[WARN] No UID found for {assembly_id}")
        return
    uid = uid_list[0]

    handle = Entrez.esummary(db="assembly", id=uid, report="full")
    record = Entrez.read(handle, validate=False)
    docsum_list = record["DocumentSummarySet"]["DocumentSummary"]
    if not docsum_list:
        print(f"[WARN] Empty DocumentSummary for {assembly_id}")
        return
    docsum = docsum_list[0]

    ftp_path = docsum.get("FtpPath_RefSeq") or docsum.get("FtpPath_GenBank")
    if not ftp_path:
        print(f"[WARN] No FTP path for {assembly_id}")
        return

    # Prefer HTTPS (more reliable than FTP)
    if ftp_path.startswith("ftp://"):
        ftp_path = ftp_path.replace("ftp://", "https://")

    base = os.path.basename(ftp_path)
    fasta_url = f"{ftp_path}/{base}_genomic.fna.gz"

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{base}_genomic.fna.gz")

    # Skip if already fully present (nonzero size)
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print(f"[SKIP] {assembly_id} already present: {out_path}")
        return

    # Download with retries; remove partial on failure
    attempts = 3
    for attempt in range(1, attempts + 1):
        try:
            print(f"[GET]  {assembly_id} -> {out_path} (try {attempt}/{attempts})")
            urllib.request.urlretrieve(fasta_url, out_path)
            if os.path.getsize(out_path) == 0:
                raise IOError("download produced zero-byte file")
            break
        except Exception as e:
            if os.path.exists(out_path):
                try:
                    os.remove(out_path)
                except Exception:
                    pass
            if attempt == attempts:
                print(f"[ERROR] {assembly_id}: {e}")
                return
            time.sleep(2)
    time.sleep(0.4)


# --- Run both projects ---
for project, fname in INPUT_FILES.items():
    print(f"\n=== Processing {project} ===")
    df = pd.read_csv(fname, sep="\t", skiprows=1, index_col=False)
    assembly_ids = df["# Assembly"].tolist()
    out_dir = os.path.join(BASE_OUTDIR, project)   # subfolder per project

    total = len(assembly_ids)
    for i, aid in enumerate(assembly_ids, 1):
        print(f"[{project} {i}/{total}] {aid}")
        try:
            fetch_fasta_from_assembly(aid, out_dir)
        except Exception as e:
            print(f"[ERROR] {aid}: {e}")
