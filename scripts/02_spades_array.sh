#SBATCH --job-name=spades_array
#SBATCH --output=/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/logs/spades_%A_%a.out
#SBATCH --error=/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/logs/spades_%A_%a.err
#SBATCH --array=1-736
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Load your conda initialization
source ~/.bashrc

# Activate the environment that contains spades.py AND pandas
micromamba activate pyenv

RAW_DIR="/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/raw_reads"
OUT_DIR="/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assembled_genomes"

cd "$RAW_DIR"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$RAW_DIR/sample_list.txt")

R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"
OUT="${OUT_DIR}/${SAMPLE}"

if [ -f "${OUT}/contigs.fasta" ]; then
    echo "[$SAMPLE] Already assembled — skipping."
    exit 0
fi

echo "[$SAMPLE] Starting SPAdes..."

spades.py \
    -1 "$R1" \
    -2 "$R2" \
    -o "$OUT" \
    --threads 4 \
    --memory 32 \
    --isolate

echo "[$SAMPLE] Finished SPAdes."
