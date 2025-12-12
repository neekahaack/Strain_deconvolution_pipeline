import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import glob
import subprocess
from Bio import AlignIO, SeqIO
import re

filenames = glob.glob("/g/typas/Personal_Folders/Neeka/Model/data/04_16s_sequences/*16S_genes.ffn")
combined_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/combined.ffn"
sequence_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/sequence_lengths.csv"
histogram_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/sequence_lengths_histogram.png"
filtered_dir = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/filtered_by_genus"
min_length = 1400
max_length = 1550
aligned_dir  = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/aligned_by_genus"
filtered_csv_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/filtered.csv"
species_count_before = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_count_before.csv"
species_count_after = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_count_after.csv"
scatter_output = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_abundance_scatter.png"
taxonomy_base = "/g/typas/Personal_Folders/Neeka/Model/data/02_taxonomy/"
matrices_dir = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/pairwise_by_genus"

#--------COMBINE ALL 16S FILES INTO 1 FILE-------------
def combine(
    filenames,
    combined_file
):
    with open (combined_file, 'w') as outfile:
    #iterate through input folder
        for names in filenames:
            with open(names) as infile:
                outfile.write(infile.read())

#------------CHECK DISTRIBUTION OF SEQUENCE LENGTHS---------------------
def counter(
    combined_file,
    sequence_file
):
    #opens file in read mode, then reads content into variable, lines, as a list of strings. 
    with open(combined_file, "r") as f:
        lines = f.readlines()

    with open(sequence_file, "w") as f:
        seq_name = ""
        seq_length = 0
        results = []
        for line in lines:
            line = line.strip() #renoves any spaces from line and rewrite back into line vcariable
        #need to iterate through input files, identify header, and count sequence
        #remember, there are multiple sequences per file sometimes
            if line.startswith(">"):
                if seq_name != "":
                    results.append((seq_name, seq_length))
                seq_name = line     # save header
                seq_length = 0      # reset length
            else:
                seq_length += len(line)

          # save the last sequence
        if seq_name != "":
            results.append((seq_name, seq_length))

        f.write("Sequence,Length\n") #column names
        for seq_name, seq_length in results:
            f.write(f"{seq_name},{seq_length}\n")

#plot distribution
def plot_distribution(
    sequence_file,
    histogram_file):
    data = pd.read_csv(sequence_file)
    sns.histplot(
        data=data, 
        x = "Length",
        bins = 30 )
    plt.tight_layout()
    plt.savefig(
        histogram_file,
        dpi=300)

#------------COMBINE TAXONOMY DOCUMENTS------------
def load_all_taxonomy_summaries(taxonomy_base):
    summary_files = glob.glob(os.path.join(taxonomy_base, "*", "gtdbtk.bac120.summary.tsv"), recursive=True)
    if not summary_files:
        raise FileNotFoundError(f"No GTDB-Tk summary files found in {taxonomy_base}")

    dfs = []
    for f in summary_files:
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
    taxonomy = pd.concat(dfs, ignore_index=True)
    return taxonomy

#----------CHECK SPECIES ABUNDANCE-----------------
def tally_species(
    taxonomy_base,
    fasta_file,
    output_csv=None,
    comparison_df=None,
    scatter_output=None
    ):
    taxonomy = load_all_taxonomy_summaries(taxonomy_base)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    ffn_df = pd.DataFrame({
        "seq_id": [rec.id for rec in records],
        "sequence": [str(rec.seq) for rec in records]})
    ffn_df["user_genome"] = ffn_df["seq_id"].str.split("_____").str[-1] #to take last element

    merged = pd.merge(ffn_df, taxonomy, on="user_genome", how="inner") #how is for rows, not columns!!

    counts = (merged["classification"].value_counts().rename_axis("classification").reset_index(name="count"))
    #counts.columns = ["classification", "Count"].  #uneccessary

    #plotting
    if comparison_df is not None and scatter_output is not None:
        merged_counts = pd.merge(
            counts,
            comparison_df,
            on="classification",
            suffixes=("_after", "_before"),
            how="outer"
        ).fillna(0)

    # Add a column showing % retained
        merged_counts["percent_retained"] = (
            merged_counts["count_after"] / merged_counts["count_before"] * 100
        )

        plt.figure(figsize=(8, 8))
        sns.scatterplot(
            data=merged_counts,
            x="count_before",
            y="count_after"
        )
        plt.title("Species abundance before vs. after filtering")
        plt.xlabel("Before filtering")
        plt.ylabel("After filtering")
        plt.xscale("log")
        plt.yscale("log")

        # Label species that changed the most (e.g. <50% retained or >200%)
        for _, row in merged_counts.iterrows():
            if row["percent_retained"] < 50 or row["percent_retained"] > 200:
                plt.text(
                    row["count_before"],
                    row["count_after"],
                    row["classification"],
                    fontsize=6,
                    alpha=0.7
                )

        plt.tight_layout()
        plt.savefig(scatter_output, dpi=300)

    # Save
    counts.to_csv(output_csv, index=False)
    return counts

#----------------FILTER SEQUENCES-------------------------
def filter_sequences(
    combined_file,
    min_length,
    max_length,
    taxonomy,
    out_dir="filtered_by_genus"
):
    # make output directory if it doesn’t exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # extract genus from classification column if not already there
    taxonomy["genus"] = taxonomy["classification"].str.extract(r"g__([^;]+)")

    # make lookup table: user_genome → genus
    genome_to_genus = taxonomy.set_index("user_genome")["genus"].to_dict()

    # read fasta sequences
    records = list(SeqIO.parse(combined_file, "fasta"))

    # iterate over records and group them by genus
    genus_records = {}  # dictionary to store list of sequences per genus
    total = 0

    for rec in records:
        total += 1
        seq_len = len(rec.seq)
        if seq_len < min_length or seq_len > max_length:
            continue  # skip if outside length range

        # get user genome from header (after last "_____")
        user_genome = rec.id.split("_____")[-1]

        # find genus
        genus = genome_to_genus.get(user_genome, None)
        if genus is None:
            continue  # skip if not in taxonomy

        # add to genus dictionary
        if genus not in genus_records:
            genus_records[genus] = []
        genus_records[genus].append(rec)

    # write one fasta file per genus
    for genus, recs in genus_records.items():
        # clean file name
        genus_filename = re.sub(r"[^A-Za-z0-9._-]", "_", genus)
        output_path = os.path.join(out_dir, f"{genus_filename}.fasta")

        with open(output_path, "w") as f:
            for rec in recs:
                f.write(f">{rec.id}\n")
                f.write(str(rec.seq) + "\n")

    print(f"Filtered {len(genus_records)} genera from {total} total sequences.")
    print(f"Files written to: {out_dir}")

    for genus, recs in sorted(genus_records.items(), key=lambda x: (-len(x[1]), x[0])):
        print(f"{genus}: {len(recs)} sequences")
        
#-------------------ALIGN ALL FILTERED GENUS FILES-----------------
def align_sequences_per_genus(
    filtered_dir,
    aligned_dir,
    threads=2
):
    # make output folder
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)

    # find input fasta files (one per genus)
    fasta_files = [f for f in os.listdir(filtered_dir) if f.endswith(".fasta")]
    print(f"Found {len(fasta_files)} genus files to align.")

    for i, file in enumerate(sorted(fasta_files)):
        input_path = os.path.join(filtered_dir, file)
        genus = file[:-6] if file.endswith(".fasta") else os.path.splitext(file)[0]
        output_path = os.path.join(aligned_dir, f"{genus}.aligned.fasta")

        # skip if already aligned and non-empty
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"[{i+1}/{len(fasta_files)}] skipping {genus} (already aligned)")
            continue

        # quick check: skip if fewer than 2 sequences
        nseq = sum(1 for _ in SeqIO.parse(input_path, "fasta"))
        if nseq < 2:
            print(f"[{i+1}/{len(fasta_files)}] skipping {genus} (nseq < 2)")
            continue

        print(f"[{i+1}/{len(fasta_files)}] aligning {genus}...")
        with open(output_path, "w") as out:
            subprocess.run([
                "mafft",
                "--thread", str(threads),
                "--auto",
                input_path
            ], stdout=out, check=True)

    print("done aligning all genera.")

#--------------COMPUTE PAIRWISE IDENTITY FOR ONE ALIGNMENT--------------
def compute_identity_matrix(
    aligned_file
):
    # read alignment
    alignment = AlignIO.read(aligned_file, "fasta")
    alignment = list(alignment)
    n = len(alignment)

    ids = [rec.id for rec in alignment]
    pairwise_matrix = pd.DataFrame(index=ids, columns=ids, dtype=float)

    # handle small cases
    if n == 0:
        return pairwise_matrix
    if n == 1:
        pairwise_matrix.iloc[0, 0] = 100.0  # self-identity
        return pairwise_matrix

    total_number_sequences = n
    for i, rec_i in enumerate(alignment):
        seq_i = str(rec_i.seq)
        pairwise_matrix.iloc[i, i] = 100.0  # diagonal = 100% identity
        print(f'Processed {i+1} out of {total_number_sequences} sequences.')
        for j in range(i + 1, n):
            rec_j = alignment[j]
            seq_j = str(rec_j.seq)

            matches = 0
            aligned = 0
            for a, b in zip(seq_i, seq_j):
                if a == "-" or b == "-":
                    continue
                aligned += 1
                if a == b:
                    matches += 1

            identity = (matches / aligned * 100.0) if aligned > 0 else 0.0
            pairwise_matrix.iloc[i, j] = identity
            pairwise_matrix.iloc[j, i] = identity

    return pairwise_matrix

#--------------BATCH: ONE MATRIX PER GENUS ALIGNMENT--------------
def compute_pairwise_per_genus(
    aligned_dir,          # folder with *.aligned.fasta (one per genus)
    matrices_dir,         # folder to write *.pairwise_identity.csv
    sep=","               # change to "\t" if you prefer TSV
):
    # make output folder
    if not os.path.exists(matrices_dir):
        os.makedirs(matrices_dir)

    aligned_files = [f for f in os.listdir(aligned_dir) if f.endswith(".aligned.fasta")]
    print(f"Found {len(aligned_files)} aligned files.")

    for i, file in enumerate(sorted(aligned_files)):
        aligned_path = os.path.join(aligned_dir, file)
        genus = file.replace(".aligned.fasta", "")
        out_path = os.path.join(matrices_dir, f"{genus}.pairwise_identity.csv")

        # skip if already computed and non-empty
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            print(f"[{i+1}/{len(aligned_files)}] skipping {genus} (matrix exists)")
            continue

        print(f"[{i+1}/{len(aligned_files)}] computing pairwise identity for {genus}...")
        try:
            matrix = compute_identity_matrix(aligned_path)
        except Exception as e:
            print(f"  -> ERROR on {genus}: {e}")
            continue

        matrix.to_csv(out_path, sep=sep)
        print(f"  -> saved: {out_path}")

    print("done computing matrices for all genera.")

if __name__ == '__main__':
    print('1) load taxonomy')
    taxonomy = load_all_taxonomy_summaries(taxonomy_base)

    print('2) filter per genus -> FASTAs in filtered_dir')
    filter_sequences(
        combined_file=combined_file,
        min_length=min_length,
        max_length=max_length,
        taxonomy=taxonomy,
        out_dir=filtered_dir
    )

    print('3) align each genus -> aligned_dir')
    align_sequences_per_genus(
        filtered_dir=filtered_dir,
        aligned_dir=aligned_dir,
        threads=2
    )

    print('4) compute pairwise identity matrices per genus -> matrices_dir')
    compute_pairwise_per_genus(
        aligned_dir=aligned_dir,
        matrices_dir=matrices_dir,
        sep=","  # use "\t" for TSV
    )

    print('All done, can I go home now?')





