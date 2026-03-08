from os.path import join
import pandas as pd
import os

configfile: "config/config.yaml"

# Read FASTQ metadata 
df = pd.read_csv(
   "config/fastqlist/fastqlist.txt",
    header="infer",
    sep="\t"
)

# Parse subject + timepoint
df["samplename"] = df["sample"]

df[["subject", "timepoint"]] = (
    df["samplename"]
    .str.split("-", expand=True)
)

# Lists used by Snakemake
samples     = sorted(df["samplename"].unique())
subjects    = sorted(df["subject"].unique())

# Define dictionary to map sample names to fastq filepaths
fastq_map = (
    df.set_index("samplename")[["fastq1", "fastq2"]]
      .to_dict(orient="index")
)

rule all:
    input:
        "workflow/out/midasOutput/species/species_profile_all.txt",
        "workflow/out/midasOutput/XAA_029/species/species_profile_subject.txt"



include: "workflow/rules/runMIDAS.smk"


