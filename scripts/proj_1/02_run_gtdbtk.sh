#!/usr/bin/env bash
# ^ shebang: tells the system to run this file using bash

#got gtdbtk program in an env: micromamba create -n gtdbtk -c conda-forge -c bioconda gtdbtk=2.* -y
#activated env: micromamba activate gtdbtk
#export GTDBTK_DATA_PATH=/scratch/nhaack/gtdbtk_db

# 0) pick a project and go there
cd /g/typas/Personal_Folders/Neeka/Model/data/assemblies/proj_number

# 1) choose ONE genome file (replace the example filename below with one from 'ls *_genomic.fna')
# ONE=/g/typas/Personal_Folders/Neeka/Model/data/assemblies/PRJNA903559/GCF_003434745.1_ASM343474v1_genomic.fna.gz

# for unzipping all files in the project
#each time through loop, variable f holds one file name with the specific ending
#[ -f "PRJNA903559_unzipped/${f%.gz}" ] checks if file in that path already exists
#|| is OR operator
for f in *.fa.gz; do
  [ -f "HiBC_Genome_sequences_20240717_unzipped/${f%.gz}" ] || \
  gunzip -c "$f" > "HiBC_Genome_sequences_20240717_unzipped/${f%.gz}"
done


# for one file: only unzip if the .fna doesn't already exist
# if [ ! -f "${ONE%.gz}" ]; then
#   gunzip -c "$ONE" > "${ONE%.gz}"
#   echo "[OK] wrote ${ONE%.gz}"
# else
#   echo "[SKIP] already have ${ONE%.gz}"
# fi

# 3) stage it in a one-file directory
# mkdir -p /g/typas/Personal_Folders/Neeka/Model/data/one_genome
# ln -sf "${ONE%.gz}" /g/typas/Personal_Folders/Neeka/Model/data/one_genome/

# 4) run GTDB-Tk on just that one file (1 CPU)tmux ls
gtdbtk classify_wf \
  --genome_dir /g/typas/Personal_Folders/Neeka/Model/data/HiBC_Genome_sequences_20240717/HiBC_Genome_sequences_20240717_unzipped \
  --out_dir    /g/typas/Personal_Folders/Neeka/Model/data/classified/HiBC_Genome_sequences_20240717 \
  --cpus 32 \
  --skip_ani_screen \
  --extension fa


set -euo pipefail
# -e : exit if any command fails
# -u : exit if we try to use an undefined variable
# -o pipefail : if we use pipelines (cmd1 | cmd2), fail if *any* part fails

# -------------------------------
# 1. Define variables (adjust these paths!)
# -------------------------------

ASSEMBLIES_ROOT="/g/typas/Personal_Folders/Neeka/Model/results/assemblies"   # where your downloaded FASTAs are
OUT_ROOT="/g/typas/Personal_Folders/Neeka/Model/results/gtdbtk_out"         # where results will go
GTDB_DATA="/scratch/nhaack/gtdbtk_db"        # where GTDB reference data should live

# -------------------------------
# 2. Make sure GTDB data exists (download if not)
# -------------------------------

if [ ! -d "$GTDB_DATA" ]; then
    echo "[INFO] GTDB data not found, downloading..."
    gtdbtk download-data --data-dir "$GTDB_DATA"
else
    echo "[INFO] GTDB data already exists at $GTDB_DATA, skipping download."
fi

# -------------------------------
# 3. Loop over all project folders in assemblies/
# -------------------------------
for PROJECT_DIR in "$ASSEMBLIES_ROOT"/*; do    #* means all subfolders in assemblies
    # skip if not a directory
    [ -d "$PROJECT_DIR" ] || continue

    PROJECT=$(basename "$PROJECT_DIR") #basename strips off path before last directory
    echo
    echo "=== Processing project: $PROJECT ==="

    # --- Uncompress FASTAs ---
    cd "$PROJECT_DIR"
    for f in *_genomic.fna.gz; do
        if [ ! -f "${f%.gz}" ]; then
            echo "[UNZIP] $f -> ${f%.gz}"
            gunzip -c "$f" > "${f%.gz}"
        else
            echo "[SKIP] already uncompressed: ${f%.gz}"
        fi
    done

    # --- Run GTDB-Tk ---
    OUTDIR="$OUT_ROOT/$PROJECT"
    mkdir -p "$OUTDIR"

    echo "[RUN] gtdbtk classify_wf on $PROJECT"
    gtdbtk classify_wf \
      --genome_dir "$PROJECT_DIR" \
      --out_dir    "$OUTDIR" \
      --data_dir   "$GTDB_DATA" \
      --cpus 8

    echo "[DONE] Results for $PROJECT -> $OUTDIR"
done