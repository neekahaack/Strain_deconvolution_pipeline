import subprocess
from pathlib import Path

#used for: simulated reads contained the correct species, abundance ratios made sense, no weird contamination

# ==========================================================
# PATHS
# ==========================================================

sim_dir = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/simulated_reads"
)

bowtie_base = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/bowtie_alignment"
)
bowtie_base.mkdir(exist_ok=True, parents=True)

metaphlan_db = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/09_samestr/db_MetaPhlAn"
)

index_name = "mpa_vJun23_CHOCOPhlAnSGB_202307"

#missed on first pass
TARGET = "Escherichia_coli_1strain_10x_ratio1_bg0"   # <-- your missed sample




# ==========================================================
# CHECK DATABASE
# ==========================================================

def ensure_metaphlan_db():
    pkl_file = metaphlan_db / f"{index_name}.pkl"
    fasta_file = metaphlan_db / f"{index_name}.fna.bz2"

    missing = [p for p in (pkl_file, fasta_file) if not p.exists()]
    if missing:
        raise RuntimeError(
            "MetaPhlAn DB incomplete.\nMissing:\n" +
            "\n".join(str(m) for m in missing)
        )

    print("MetaPhlAn DB found:")
    print(" ", pkl_file)
    print(" ", fasta_file)


# ==========================================================
# RUN METAPHLAN
# ==========================================================

def run_metaphlan_for_all(threads=8):
    ensure_metaphlan_db()

    

    for subdir in sorted(sim_dir.iterdir()):
        if not subdir.is_dir():
            continue

         #for missing sample only
        if subdir.name != TARGET:
            continue

        metagenome_name = subdir.name
        print(f"\n=== Running MetaPhlAn for {metagenome_name} ===")

        # find FASTQs
        r1_list = sorted(subdir.glob("*_1.fq"))
        if not r1_list:
            print(f"No R1 FASTQ found in {subdir}, skipping.")
            continue

        r1 = r1_list[0]
        r2 = subdir / r1.name.replace("_1.fq", "_2.fq")
        if not r2.exists():
            print(f"Missing R2 for {metagenome_name}, skipping.")
            continue

        # output paths
        out_dir = bowtie_base / metagenome_name
        out_dir.mkdir(exist_ok=True, parents=True)

        profile_out = out_dir / f"{metagenome_name}.profile.txt"
        #sam_out = out_dir / f"{metagenome_name}.sam.bz2"
        bowtie2out = out_dir / f"{metagenome_name}.metaphlan.bowtie2.bz2"

        #if profile_out.exists():
        #    print(f"Already exists: {profile_out}")
        #    continue

        # proper paired-end input
        fastq_input = f"{r1},{r2}"

        cmd = [
            "metaphlan",
            fastq_input,
            "--input_type", "fastq",
            "--db_dir", str(metaphlan_db),
            "-x", index_name,
            "--mapout", str(bowtie2out),
            #"-s", str(sam_out),
            "--offline",
            "--nproc", str(threads),
            "-t", "rel_ab",
            "-o", str(profile_out),
        ]

        print("COMMAND:\n ", " ".join(cmd))
        subprocess.run(cmd, check=True)

        print(f"MetaPhlAn profile created: {profile_out}")
        #print(f"MetaPhlAn SAM created:     {sam_out}")


# ==========================================================
# MAIN
# ==========================================================

if __name__ == "__main__":
    run_metaphlan_for_all(threads=16)
    print("\nAll MetaPhlAn profiling complete!\n")


