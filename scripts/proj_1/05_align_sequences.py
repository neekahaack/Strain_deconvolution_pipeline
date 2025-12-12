import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import glob
import subprocess
from Bio import AlignIO, SeqIO

filenames = glob.glob("/g/typas/Personal_Folders/Neeka/Model/data/04_16s_sequences/*16S_genes.ffn")
combined_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/combined.ffn"
sequence_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/sequence_lengths.csv"
histogram_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/sequence_lengths_histogram.png"
filtered_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/filtered.ffn"
min_length = 1400
max_length = 1550
aligned_file  = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/aligned_mafft_auto.ffn"
#filtered_csv_file = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/filtered.csv"
species_count_before = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_count_before.csv"
species_count_after = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_count_after.csv"
scatter_output = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/species_abundance_scatter.png"
taxonomy_base = "/g/typas/Personal_Folders/Neeka/Model/data/02_taxonomy/"
output_matrix_path = "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/pairwise_identity_16s.csv"

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

#------------CHECK DISTRIBUTION---------------------
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
    #Find and combine all gtdbtk.bac120.summary.tsv files in subfolders
    #one star tells it to not go nore than one layer deep with search
    summary_files = glob.glob(os.path.join(taxonomy_base, "*", "gtdbtk.bac120.summary.tsv"), recursive=True)
    if not summary_files:
        raise FileNotFoundError(f"No GTDB-Tk summary files found in {taxonomy_base}")

    dfs = []
    for f in summary_files:
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
    taxonomy = pd.concat(dfs, ignore_index=True)
    return taxonomy

#to check
taxonomy = load_all_taxonomy_summaries(taxonomy_base)
len(taxonomy)

#----------CHECK SPECIES ABUNDANCE-----------------
#count species abundance
#this is prior to filtering
#but also want to be able to call function after filtering on filtered_file
def tally_species(
    #species_count,
    taxonomy_base,
    fasta_file,
    #combined_file
    output_csv=None,
    comparison_df=None,
    scatter_output=None
    ):
    #load taxonomy data into variable taxonomy
    taxonomy = load_all_taxonomy_summaries(taxonomy_base)
    #parse file with sequences
    records = list(SeqIO.parse(fasta_file, "fasta"))
    # Extract header info
    ffn_df = pd.DataFrame({
        "seq_id": [rec.id for rec in records],
        "sequence": [str(rec.seq) for rec in records]})
    #put dataframe into .csv
    ffn_df["user_genome"] = ffn_df["seq_id"].str.split("_____").str[-1] #to take last element
    #ffn_df.to_csv(filtered_csv_file, index = False)       #not necessary to convert file type before merging
    #print("taxonomy columns:", taxonomy.columns.tolist())
    #print("ffn_df columns:", ffn_df.columns.tolist())
    
    # Merge #on=column they have in common #
    merged = pd.merge(ffn_df, taxonomy, on="user_genome", how="inner") #how is for rows, not columns!!

    # Tally
    counts = merged["classification"].value_counts().reset_index() #not sure if the reset_index is necessary
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
#need to do some kind of iteration that goes through, detects header, then goes to next line
#count characters in next line, header and next line into new file if character count exceeds certain number

def filter_sequences(
    combined_file,
    filtered_file,
    min_length,
    max_length
):
    with open(combined_file, "r") as f:
        lines = f.readlines()

    with open(filtered_file, "w") as f:
        header = ""
        seq = ""
        for line in lines:
            line = line.strip() #renoves any spaces from line and rewrite back into line variable
        #need to iterate through input lines, identify header, and count sequence
        #remember, there are multiple sequences per file sometimes
            if line.startswith(">"):
                if len(seq) >= min_length and len(seq) <= max_length:
                    f.write(header + "\n")
                    f.write(seq + "\n")
                header = line
                seq = ""
            else:
                seq += line

          # save the last sequence
        if len(seq) >= min_length and len(seq) <=max_length:
            f.write(header + "\n")
            f.write(seq + "\n")

#-------------------ALIGN THE SEQUENCES-----------------
#problem is that mafft has no python wrapper or module like prokka or gtdbtk, so have to run with command line
#or do a subprocess
#then per sequence completeness score in same file

# run mafft and save alignment
def align_sequences(
    filtered_file,
    aligned_file
):
    with open(aligned_file, "w") as out:    # temporary variable to refer to file object (NOT the files contents)
        subprocess.run([                    #with is context manager, so we open and it closes the file for us
            "mafft", 
            "--thread", 
            "2",            #cpu
            "--auto",          # let mafft choose the method, without auto its doing default = FFT-NS-2, which auto is doing too bc of dataset size
            filtered_file        # being run on the filter file
        ], stdout=out)             #write the standard output to temporary variable



#--------------COMPUTE PAIRWISE IDENTITY----------------
#build empty matrix, need to know how big to make it so need to see how many sequences make it into aligned file
# read the alignment is a way
def compute_identity(
    aligned_file
):
    alignment = AlignIO.read(aligned_file, "fasta")  #to variable alignment
    alignment = list(alignment)
    # Below jus for debugging
    #alignment = alignment[0:25]
#n = len(alignment)
#print(n)  #gives 3081 when min_length is 1400, need to rerun script to see what happens when max is set
#but maybe this isnt necessary, because the alignment already has rec.id and rec.seq, so it aready knows the length
#should make empty 3081 x 3081 matrix with headers 
# use pairwise_matrix = pd.DataFrame
    pairwise_matrix = pd.DataFrame(
        index = [rec.id for rec in alignment],
        columns = [rec.id for rec in alignment])

#assign every sequence in column and row headers
#loop, nic said twice but not sure why, maybe first is assigning location and second is doing the calculation?
#so that every index is pairwise allignment score for each pair
#then not sure how to extract important values from that but thats a problem for later
    total_number_sequences = len(alignment)

    for i, rec_i in enumerate(alignment):   #enumerate to go through indices and element at each
        seq_i = str(rec_i.seq)              #convert seq to string
        print(f'Processed {i} out of {total_number_sequences} sequences.')
        for j, rec_j in enumerate(alignment):
            if j <= i: #this is so only upper triangle above diagonal gets computed
                continue  #means exit loop at this point and move onto next loop
            seq_j = str(rec_j.seq)          #convert seq to string
            matches, aligned = 0, 0         #variables values initially set to 0 and 0
            for a, b in zip(seq_i, seq_j):     #zip is python function, zip pairs each base of seq_i and seq_j position by position
                if a == "-" or b == "-":        #when looping, if a and b are gaps
                    continue        #means exit loop at this point and move onto next iteration, aka skip this position and go to next
                aligned += 1                    #ELSE, add 1 to aligned counter
                if a == b:                      #if sequence matches
                    matches += 1                #add 1 to matches counter
            identity = (matches / aligned * 100) if aligned > 0 else 0
            pairwise_matrix.iloc[i, j] = identity
            pairwise_matrix.iloc[j, i] = identity       #write into both haves of matrix, even tho both haves are identitical.
    return pairwise_matrix
# First we loop through the alignment file. 
# At every index (i), we take one sequence record (rec_i) and convert its sequence to a string.
# For every iteration through that outer loop, we also do another full iteration through the same alignment (the inner loop).
# In that inner loop, we skip the index if it's less than or equal to the outer loop index (j <= i).
# This ensures we only compute and populate the matrix for the upper triangle (no self-comparisons or duplicates).
# Inside those two loops, we compare the two sequences (seq_i and seq_j) position by position using the zip() function.
# For each pair of characters (a, b) from seq_i and seq_j:
#     - If either a or b is a gap ("-"), we skip this position (continue to next).
#     - Otherwise, we add +1 to the aligned counter.
#     - And if a and b are exactly the same (a == b), we also add +1 to the matches counter.
# Once the loop over positions is done, we calculate the percent identity as (matches / aligned * 100).
# This calculation happens inside the inner loop (for each pair of sequences).
# Then we store that identity value into the DataFrame (pairwise_matrix)
# at the position [i, j] and [j, i], using pairwise_matrix.iloc[i, j] = identity.
#i in .iloc is for integer. 

#then to figure out what those nonconsensus species are
#load taxonomy file, group by species, pipe into tally() to count how many per species, before and after filtering for 16S length (try kicking out over 1550), then scatter the counts
#this is scatterplot before and after, use full.join. 



# # rewrite with identity percentages
#     with open(aligned_file, "w") as out:                #open file, it gets EMPTIED, and assign it to handle out
#         for record in alignment:                         #each record has .id and .seq from previous alignment variable!
#             seq = str(record.seq)                       #take .seq portion and convert to string
#             non_gaps = sum(1 for base in seq if base != "-")    #count all things that arent gaps 
#             total = len(seq)                                    #count total sequence length
#             identity = (non_gaps / total) * 100 if total > 0 else 0     #do the calculation
#             out.write(f">{record.id}\n")                               #rewrite the .id
#             out.write(seq + "\n")              #rewrite the seq, could also be written as out.write(f"{record.seq}\n") 
#             out.write(f"Identity: {identity:.2f}%\n")       #write percentage with 2 decimal points    

#jalview can tell me if there is one sequence driving gaps

if __name__ == '__main__':
    # print('1')
    # combine(
    #     filenames,
    #     combined_file)
    # print('2')
    # counter(
    #     combined_file,
    #     sequence_file)
    # print('3')
    # plot_distribution(
    #     sequence_file,
    #     histogram_file)
    # print('4')
    # before_counts = tally_species(
    #     taxonomy_base=taxonomy_base,
    #     fasta_file=combined_file,
    #     output_csv=species_count_before)
    # print('5')
    # filter_sequences(
    #     combined_file,
    #     filtered_file,
    #     min_length,
    #     max_length)
    # print('6')
    # #basically want to check the species abundance before and after filtering
    # after_counts = tally_species(
    #     taxonomy_base=taxonomy_base,
    #     fasta_file=filtered_file,
    #     output_csv=species_count_after,
    #     comparison_df=before_counts,
    #     scatter_output=scatter_output)
    #align_sequences(
    #    filtered_file,
    #    aligned_file)
    output = compute_identity(
        aligned_file)
    # TODO: Write output to disk
    output.to_csv(output_matrix_path)
    print('All done, can I go home now?')