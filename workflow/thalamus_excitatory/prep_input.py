import os
import sys
import subprocess
import re
import numpy as np
import pandas as pd

# Specify a path to input metadata table
meta_table = "../../input/thalamus_excitatory/combined_thalamus_metadata.csv"
# Specify a path to input bam files
bam = 'possorted_genome_bam.bam'
# Specify output directory
outdir = "../../input/thalamus_excitatory/bam"
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Open the metadata
df = pd.read_csv(meta_table, index_col="SampleID")

for id in df.index:
    print(id, "being processed...")
    # Specify paths to bam and bam.bai
    prefix = df.loc[id, 'path']
    files = [f"{prefix}/{id}/{bam}", f"{prefix}/{id}/{bam}.bai"]

    for file in files:
        if not os.path.exists(file):
            raise ValueError(f"{file} not found!")
        else:
            disease_study = df.loc[id, 'study_specific_disease_specific']
            newid = f"{id}_{disease_study}"
            filename = file.split("/")[-1]
            # Specify command
            outfile = f"{outdir}/{newid}_{filename}"
            cmd = f"cp {file} {outfile}"
            # Copy
            if not os.path.exists(outfile):
                print(file, "copied to", outfile)
                subprocess.run(cmd, shell=True)
            else:
                print(file, "already copied. skipped!")




# print(df)

