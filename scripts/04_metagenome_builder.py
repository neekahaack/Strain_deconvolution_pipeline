import pandas as pd
from pathlib import Path
import json
import itertools
import random

# ==========================================================
#                 PARAMETER SPACE
# ==========================================================

PARAMS = {
    #Species that you want to simulate as "strain(s) of interest"
    "species_of_interest": [
        #"Lachnospiraceae nov.",
        "Escherichia coli"
        #'Bacteroides uniformis'
    ],

    #How many conspecific strains to simulate together (1, 2, 3)
    "conspecific_counts": [1, 2, 3],

    #Main strain coverages (depth of sequencing)
    "target_coverages": [
        1, 
        2, 
        3,
        4,
        5, 
        6, 
        7,
        8,
        9,
        10, 
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20, 
        100, 
        200],

    #Ratios of secondary strains relative to the main one
    "ratios": [1, 0.5, 0.2, 0.1, 0.01, 0.001, 0.0001],

    #Background microbial complexities
    "background_complexities": [0, 2, 5, 10, 20, 50],

    #Coverage for background microbes
    "background_coverage": 20,

    #Random seed for reproducibility
    "random_seed": 42,
}

# ==========================================================
#                   PATHS
# ==========================================================

metadata_tsv = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/assemblies_metadata.tsv"
)

output_config_dir = Path(
    "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/metagenome_configs"
)
output_config_dir.mkdir(exist_ok=True)


# ==========================================================
#                 HELPER FUNCTIONS
# ==========================================================

def pick_random(df, n):
    """Pick n random rows."""
    return df.sample(n=n)


def build_config_name(sp, n_strains, cov, r, bg):
    """Make readable config filenames."""
    return f"{sp.replace(' ', '_')}_{n_strains}strain_{cov}x_ratio{r}_bg{bg}"


def build_config(sample_id, genomes):
    """Builds a simulation config dictionary."""
    return {
        "sample_id": sample_id,
        "genomes": genomes
    }


# ==========================================================
#               LOAD METADATA
# ==========================================================

df = pd.read_csv(metadata_tsv, sep="\t")
#Initialize coverage column (needed for assignment later)
df["coverage"] = 0

random.seed(PARAMS["random_seed"])

configs = []

# ==========================================================
#            MAIN PARAMETER SPACE LOOPS
# ==========================================================

for species in PARAMS["species_of_interest"]:

    df_sp = df[df["Species"] == species]   #keeps rows where value=true is used to filter, uses boolean series

    if df_sp.empty:
        print(f"⚠ No isolates found for species: {species}. Skipping.")
        continue

    for n_strains in PARAMS["conspecific_counts"]:
        if len(df_sp) < n_strains:
            print(f"⚠ Not enough conspecific strains for {Species} (need {n_strains}). Skipping.")
            continue

        for cov in PARAMS["target_coverages"]:
            for r in PARAMS["ratios"]:

                #Select conspecific strains
                #####################################################################
                # TODO: This is a quick and dirty fix - make this proper moving ahead
                #####################################################################
                
                #conspecific_group = df_sp.sample(n_strains, random_state=PARAMS["random_seed"]).copy()
                conspecific_group = df_sp.iloc[:n_strains, :]

                #Coverage assignment
                #Assign coverage to conspecific strains
                conspecific_group.loc[conspecific_group.index[0], "coverage"] = cov
                conspecific_group.loc[conspecific_group.index[1:], "coverage"] = cov * r

                #Prepare background pool
                df_bg_pool = df[~df["Sample_id"].isin(conspecific_group["Sample_id"])]

                for bg in PARAMS["background_complexities"]:

                    #Choose background species
                    if bg > 0:
                        if len(df_bg_pool) < bg:
                            print(f"⚠ Not enough genomes for background complexity {bg}. Skipping.")
                            continue
                        bg_genomes = pick_random(df_bg_pool, bg).copy()
                        bg_genomes["coverage"] = PARAMS["background_coverage"]
                    else:
                        bg_genomes = pd.DataFrame(columns=df.columns)

                    #Combine all genomes
                    combined = pd.concat([conspecific_group, bg_genomes], ignore_index=True)

                    #Convert to list of dicts
                    genome_entries = [
                        {
                            "id": row["Sample_id"],
                            "species": row["Species"],
                            "fasta": row["Assembly_path"],
                            "coverage": row["coverage"]
                        }
                        for _, row in combined.iterrows()
                    ]

                    #Create config
                    name = build_config_name(species, n_strains, cov, r, bg)
                    configs.append(build_config(name, genome_entries))


# ==========================================================
#                   SAVE CONFIGS
# ==========================================================

for cfg in configs:
    out_path = output_config_dir / f"{cfg['sample_id']}.json"
    with open(out_path, "w") as f:
        json.dump(cfg, f, indent=2)

print(f"Generated {len(configs)} metagenome configs!")
