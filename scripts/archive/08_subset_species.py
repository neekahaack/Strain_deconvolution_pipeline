import subprocess
from pathlib import Path        #use Path instead of just putting file pathway in as a string
import pandas as pd
import requests

#--------------CHOOSE SUBSET OF SPECIES FROM SUPPLEMENTARY TABLE 1----------
#choose 10 species
df = pd.read_csv(
    "/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/41587_2018_9_MOESM3_ESM.csv",
    sep=";",
    skiprows = 1)

#keep only HBC isolates (the ones cultured and sequenced by the study), the other ones are from public datasets
df = df[df["Source"] == "HBC"]

#subset 10 randomly with only 1 of each species
subset = df.drop_duplicates(subset="Species").sample(n=10, random_state=42)     #this is called a seed (the random state thing)

#save to a new CSV for reference
subset.to_csv("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/selected_10_species.csv", index=False)

print(subset[["Identifier", "Genus", "Species"]])


