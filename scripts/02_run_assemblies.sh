set -e

RAW_DIR="/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/raw_reads"
ASM_DIR="/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assembled_genomes"
CSV_PATH="/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/41587_2018_9_MOESM3_ESM.csv"
META_OUT="/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assemblies_metadata.tsv"

echo "Ensuring assembled_genomes directory exists..."
mkdir -p "$ASM_DIR"

cd "$RAW_DIR"

echo "Creating sample_list.txt ..."
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > sample_list.txt
NUM=$(wc -l < sample_list.txt)
echo "Found $NUM samples."

echo "Making log directory ..."
mkdir -p /g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/logs

echo "Submitting SPAdes job array ..."
JOBID=$(sbatch --array=1-$NUM /g/typas/Personal_Folders/Neeka/Model/scripts/spades_array.sh | awk '{print $4}')
echo "Submitted as job ID: $JOBID"

echo "Waiting for SPAdes array to finish ..."
while squeue -j "$JOBID" 2>/dev/null | grep -q "$JOBID"; do
    sleep 20
done

echo "All SPAdes jobs finished."

