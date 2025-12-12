#first thing to look at could be to check if theres grouping within genera
#so headers in matrix are gene id___seq id
#need to load taxonomy file, maybe can use from last script since I did "return taxonomy"? 
#then with classification column to separate species, genus, etc
#these can be kind of linked (?) 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

taxonomy_base = "/g/typas/Personal_Folders/Neeka/Model/data/02_taxonomy/"
pivoted_matrix = "/g/typas/Personal_Folders/Neeka/Model/data/06_analysis/pairwise_identity_16s_pivoted_long.csv"


#-------------IMPORT NECESSARY FILES---------------
#import combined taxonomy
align_sequences = __import__('05_align_sequences')
taxonomy = align_sequences.load_all_taxonomy_summaries(taxonomy_base)
len(taxonomy)
#can check it by writing in ipython terminal : len(texonomy)

#import matrix
matrix = pd.read_csv('/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/pairwise_identity_16s.csv', index_col=0)
#can check it by writing in ipython terminal : matrix.iloc[:10, :10]
#there will be NaN now instead of just a space bc pandas reads empty cells like that

#extract genus from classification 
taxonomy["genus"] = taxonomy["classification"].str.extract("g__([^;]+)")
taxonomy.head()
#so at this point, taxonomy has user_genome column with just genome id
#but pairwise identity matrix has has seq_id____user_genome (because there were multiple 16s genes for some genomes)
#make sure matrix and taxonomy have overlapping IDs, this makes A NEW COLUMN WITH JUST USER_GENOME, doesnt change headers
matrix["user_genome"] = matrix.index.str.split("_____").str[-1]

#build a lookup table: genome → genus
genome_to_genus = taxonomy.set_index("user_genome")["genus"]

#map that lookup onto the matrix rows, makes new genus column on matrix
matrix["genus"] = matrix["user_genome"].map(genome_to_genus)

#see if mapping worked, both columns exist
matrix[["user_genome", "genus"]].head()
#if 0, means every sequence has valid genus
matrix["genus"].isna().sum()

#-----------------RESHAPE DATAFRAMES--------------
#Need to put matrix in long form, do this melt thing
long_df = (
    matrix
    .drop(columns=["user_genome", "genus"])       #remove columns or rows
    .reset_index(names="seq_i")                  #turn row labels back into a column
    .melt(id_vars=["seq_i"], var_name="seq_j", value_name="identity") #wide→long format
)
long_df.head()

# Attach genus info for both sides
long_df["genus_i"] = long_df["seq_i"].map(matrix["genus"]) #replace values based on lookup
long_df["genus_j"] = long_df["seq_j"].map(matrix["genus"])

# Drop self-comparisons and NaNs
long_df = long_df.dropna(subset=["identity"])
long_df = long_df[long_df["seq_i"] != long_df["seq_j"]]

long_df.to_csv(pivoted_matrix)
#to check long_df.columns, gets Index(['seq_i', 'seq_j', 'identity', 'genus_i', 'genus_j'], dtype='object')

#-----------------INTRAGENUS COMPARISON--------------
#keep only pairs where sequences are from same genus
intragenus = long_df[long_df["genus_i"] == long_df["genus_j"]]

#check stats
intragenus_stats = (
    intragenus
    .groupby("genus_i")["identity"]
    .agg(["mean", "median", "count"])   #aggrete to compute summary stats
    .reset_index()  #reshape results into clean summary table
    .rename(columns={"genus_i": "genus"})
)
intragenus_stats.head()

#-------------PLOTTING-------------
#scatterplot
sns.scatterplot(
    data=intragenus_stats,
    x="count",
    y="median",
    size="count",
    hue="median",
    palette="viridis",
    alpha=0.8,
    legend=False
)
plt.xlabel("Number of pairwise comparisons (count)")
plt.ylabel("Median pairwise identity (%)")
plt.title("Intragenus 16S rRNA identity across genera")
plt.tight_layout()
plt.savefig(
    "/g/typas/Personal_Folders/Neeka/Model/data/06_analysis/intragenus_identity_scatter.png",
    dpi=300,
    bbox_inches="tight"
)

#bar plot
top = intragenus_stats.sort_values("median", ascending=False)
plt.figure(figsize=(8, 16))
sns.barplot(
    data=top,
    y="genus",
    x="median",
    palette="crest"
)
plt.xlabel("Median pairwise identity (%)")
plt.ylabel("Genus")
plt.title("Intragenus 16S rRNA sequence identity")
plt.tight_layout()
plt.savefig(
    "/g/typas/Personal_Folders/Neeka/Model/data/06_analysis/intragenus_identity_bar.png",
    dpi=300,
    bbox_inches="tight"
)