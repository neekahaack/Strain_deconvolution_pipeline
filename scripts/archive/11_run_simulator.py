import subprocess
from pathlib import Path
import pandas as pd
import requests

# ---Simulate reads from each genome---
genome_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assembled_genomes")
output_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/simulated_reads")
output_dir.mkdir(exist_ok=True)

#Find all contigs.fasta files in subdirectories
for genome in genome_dir.glob("*/contigs.fasta"):           #* one star makes it 1 level deep for searching
    base = genome.parent.name               # e.g. ERR2221104
    out_prefix = output_dir / f"{base}_"

    print(f"Simulating reads for {base}...")

    cmd = [
        "art_illumina",
        "-ss", "HS25",      # HiSeq 2500
        "-i", str(genome),
        "-p",               # paired-end
        "-l", "100",        # read length
        "-f", "100",        # coverage
        "-m", "500",        # mean insert size
        "-s", "50",         # std dev
        "-o", str(out_prefix)
    ]
    subprocess.run(cmd, check=True)
    print(f"finished base")

print("Simulation complete!")

# --- NOTES ---
# Illumina model and parameters (HiSeq 2500, 2×100 bp, ~500±50 bp insert) derived from ENA metadata.
# For actual benchmarking, ensure your genomes and reads correspond to the same isolate source.
 #output file 1 and 2 are the 2 directions of pair, and aln says where simulation was modeled off of

    #chose flags based on this [haack@login2 Neeka]$ curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR2221356&result=read_experiment&fields=library_name,library_layout,library_strategy,library_source,library_selection,instrument_model,nominal_length,nominal_sdev,read_count,base_count&format=tsv"
#experiment_accession    library_layout  library_strategy        library_source  library_selection       instrument_model        nominal_length  nominal_sdev   read_count      base_count      library_name    run_accession
#ERX2274729      PAIRED  WGS     GENOMIC RANDOM  Illumina HiSeq 2500     500             5298350 662293750       DN447705L:F3    ERR2221356

#cat ERR1203919_1.fastq | head -400 | awk '(NR%4==2){print length($0)}' | sort | uniq -c
#    100 100                    means 100 reads in first 400 lines, read length is 100 bp