import subprocess
from pathlib import Path
import gzip
import shutil

# ==========================================================
# PATHS
# ==========================================================

samestr_db = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/09_samestr/samestr_db"
)

db_markers = samestr_db / "db_markers"

marker_fasta = samestr_db / "samestr_markers.fna"

index_prefix = samestr_db / "samestr_marker_index"

sim_dir = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/simulated_reads"
)

out_base = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/bowtie_alignment"
)
out_base.mkdir(exist_ok=True)


# ==========================================================
# 1) BUILD SameStr MARKER FASTA
# ==========================================================

def build_samestr_marker_fasta():
    if marker_fasta.exists():
        print(f"Marker FASTA already exists: {marker_fasta}")
        return

    print(f"Building combined SameStr marker FASTA at: {marker_fasta}")

    with open(marker_fasta, "wb") as out_f:
        for gz_fa in sorted(db_markers.rglob("*.fa.gz")):
            print(f"  adding {gz_fa}")
            with gzip.open(gz_fa, "rb") as f_in:
                shutil.copyfileobj(f_in, out_f)

    print("Done building marker FASTA.")


# ==========================================================
# 2) BUILD Bowtie2 INDEX
# ==========================================================

def build_bowtie_index():
    # Bowtie2 LARGE indexes always produce *.1.bt2l
    bt2l_file = Path(f"{index_prefix}.1.bt2l")

    if bt2l_file.exists():
        print(f"Bowtie2 index already exists at: {index_prefix}")
        return

    print("Building Bowtie2 index for SameStr markers...")
    subprocess.run(["bowtie2-build", str(marker_fasta), str(index_prefix)], check=True)
    print("Index built!")

# ==========================================================
# 3) ALIGN SIMULATED READS
# ==========================================================

def align_all_simulated_reads():

    for subdir in sorted(sim_dir.iterdir()):
        if not subdir.is_dir():
            continue

        metagenome_name = subdir.name
        print(f"\n=== Processing metagenome: {metagenome_name} ===")

        out_dir = out_base / metagenome_name
        out_dir.mkdir(exist_ok=True)

        # Find read pairs
        for r1 in sorted(subdir.glob("*_1.fq")):
            base = r1.name.replace("_1.fq", "")
            r2 = subdir / f"{base}_2.fq"

            if not r2.exists():
                print(f"Missing R2 for {base}, skipping.")
                continue

            bam_sorted = out_dir / f"{base}.bam"

            if bam_sorted.exists():
                print(f"Skipping {base}: BAM already exists.")
                continue
 
            print(f"  Aligning {base}")

            sam = out_dir / f"{base}.sam"
            bam_unsorted = out_dir / f"{base}.unsorted.bam"
            

            # bowtie2
            subprocess.run([
                "bowtie2",
                "-x", str(index_prefix),
                "-1", str(r1),
                "-2", str(r2),
                "--threads", "8",
                "-S", str(sam)
            ], check=True)

            # SAM → BAM
            subprocess.run(["samtools", "view", "-bS", "-o", str(bam_unsorted), str(sam)], check=True)

            # Sort BAM
            subprocess.run(["samtools", "sort", "-o", str(bam_sorted), str(bam_unsorted)], check=True)

            # Index BAM
            subprocess.run(["samtools", "index", str(bam_sorted)], check=True)

            # Cleanup
            sam.unlink()
            bam_unsorted.unlink()

            print(f" Finished {base}")


# ==========================================================
# MAIN
# ==========================================================

if __name__ == "__main__":
    build_samestr_marker_fasta()
    build_bowtie_index()
    align_all_simulated_reads()
    print("\nAll alignments complete.")
