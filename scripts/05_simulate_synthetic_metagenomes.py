import json
import subprocess
from pathlib import Path
import shutil

# ==========================================================
#   SELECT WHICH METAGENOME CONFIGS TO SIMULATE
#   (EDIT THIS LIST MANUALLY)
# ==========================================================

SELECTED_CONFIGS = [
    #"Escherichia_coli_1strain_10x_ratio1_bg0.json",
    #"Escherichia_coli_2strain_20x_ratio0.5_bg5.json",
    #"Lachnospiraceae_nov._3strain_50x_ratio0.1_bg10.json",
    #"Lachnospiraceae_nov._2strain_100x_ratio1_bg0.json"

    ### Question 1: Assuming the simplest scenario (no bg, 1 strain, ...), what is minimum coverage required for sameStr to be able to call strain identity
    # "Escherichia_coli_1strain_200x_ratio1_bg0.json", # this would give us ground truth SNP profile
     "Escherichia_coli_1strain_100x_ratio1_bg0.json", # this would also give us ground truth SNP profile
    # "Escherichia_coli_1strain_20x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_19x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_18x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_17x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_16x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_15x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_14x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_13x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_12x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_11x_ratio1_bg0.json",
     "Escherichia_coli_1strain_10x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_9x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_8x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_7x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_6x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_5x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_4x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_3x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_2x_ratio1_bg0.json",
    # "Escherichia_coli_1strain_1x_ratio1_bg0.json",
    # ### Question 2: Same as above, but with 2 or 3 strains of the same species present at the same time
    # "Escherichia_coli_2strain_20x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_20x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_19x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_19x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_18x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_18x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_17x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_17x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_16x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_16x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_15x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_15x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_14x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_14x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_13x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_13x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_12x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_12x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_11x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_11x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_10x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_10x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_9x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_9x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_8x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_8x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_7x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_7x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_6x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_6x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_5x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_5x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_4x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_4x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_3x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_3x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_2x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_2x_ratio1_bg0.json",
    # "Escherichia_coli_2strain_1x_ratio1_bg0.json",
    # "Escherichia_coli_3strain_1x_ratio1_bg0.json"
    #...

]

CONFIG_DIR = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/metagenome_configs")
OUTPUT_BASE = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/simulated_reads")
OUTPUT_BASE.mkdir(exist_ok=True)

# ==========================================================
#                 HELPER: RUN COMMAND
# ==========================================================

def run(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

# ==========================================================
#              SIMULATE ONE METAGENOME
# ==========================================================

def simulate_metagenome(config_path, output_dir, read_length=100, insert_mean=500, insert_sd=50, print_command_only = False):

    # Load JSON config
    with open(config_path, "r") as f:
        config = json.load(f)

    sample_id = config["sample_id"]
    genomes = config["genomes"]

    if not print_command_only:
        print(f"\n=== Simulating metagenome: {sample_id} ===")
        print(f"Genomes included: {len(genomes)}")

    # Create per-metagenome temp directory
    temp_dir = output_dir / f"{sample_id}_temp"
    temp_dir.mkdir(exist_ok=True)

    # Final merged fastqs
    out_r1 = output_dir / f"{sample_id}_1.fq"
    out_r2 = output_dir / f"{sample_id}_2.fq"

    # Overwrite old outputs
    if out_r1.exists(): out_r1.unlink()
    if out_r2.exists(): out_r2.unlink()

    # -----------------------------------------------------
    # Simulate reads for each genome (ART)
    # -----------------------------------------------------
    for genome in genomes:
        genome_id = genome["id"]
        fasta = genome["fasta"]
        cov = genome["coverage"]

        prefix = temp_dir / genome_id

        if not print_command_only:
            print(f"\nSimulating {genome_id} at {cov}x coverage...")

        cmd = [
            "art_illumina",
            "-ss", "HS25",
            "-i", fasta,
            "-p",                   # paired-end
            "-l", str(read_length),
            "-f", str(cov),
            "-m", str(insert_mean),
            "-s", str(insert_sd),
            "-o", str(prefix),
        ]

        if (print_command_only):
            print(' '.join(cmd))
            return

        else:
            run(cmd)

    # -----------------------------------------------------
    # Merge all genome reads into final FASTQs
    # -----------------------------------------------------
    print("\nMerging reads...")

    with open(out_r1, "wb") as fout1, open(out_r2, "wb") as fout2:
        for genome in genomes:
            prefix = temp_dir / genome["id"]

            r1 = Path(str(prefix) + "1.fq")
            r2 = Path(str(prefix) + "2.fq")

            if r1.exists():
                fout1.write(r1.read_bytes())
            if r2.exists():
                fout2.write(r2.read_bytes())

    print(f"Final metagenome FASTQs:\n  {out_r1}\n  {out_r2}")

    # -----------------------------------------------------
    # Cleanup
    # -----------------------------------------------------
    shutil.rmtree(temp_dir)
    print(f"Removed temp directory: {temp_dir}")

    print(f"\nCompleted: {sample_id}\n")


# ==========================================================
#               MAIN: RUN SELECTED CONFIGS
# ==========================================================

if __name__ == "__main__":

    for cfg_name in SELECTED_CONFIGS:
        print(cfg_name)
        cfg_path = CONFIG_DIR / cfg_name

        if not cfg_path.exists():
            print(f"⚠ Config not found: {cfg_path}")
            continue

        # Make a dedicated output folder for this metagenome
        out_dir = OUTPUT_BASE / cfg_path.stem
        out_dir.mkdir(exist_ok=True)

        simulate_metagenome(cfg_path, out_dir, print_command_only = False) #changed to false so it actually runs art

    print("\nAll selected simulations finished!")


#Per metagenome:
# Run ART for each genome
# ↓
# Get per-genome FASTQs
# ↓
# Merge all genome-specific FASTQs into one R1 and one R2
# ↓
# Use those two files for Bowtie2, BWA, downstream steps