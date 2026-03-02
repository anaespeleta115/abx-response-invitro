from os.path import join
import pandas as pd
import os


# Read FASTQ metadata table
# ----------------------------
df = pd.read_csv(
   "config/fastqlist/fastqlist.txt",
    header="infer",
    sep="\t"
)

# sanity check (VERY helpful)
print(df.head())


df = df.rename(columns={
    "read1": "fastq1",
    "read2": "fastq2"
})


# Parse subject + timepoint
# ----------------------------
df["samplename"] = df["sample"]

df[["subject", "timepoint"]] = (
    df["samplename"]
    .str.split("-", expand=True)
)


# Lists used by Snakemake
# ----------------------------
samples     = sorted(df["samplename"].unique())
subjects    = sorted(df["subject"].unique())


rule all:
    input:
        #expand("/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/HouseholdTransmission-Stool-{sample}-trimmed-pair1.fastq.gz",sample=samplelanes),
        #expand("/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/HouseholdTransmission-Stool-{sample}.done",sample=samplelanes)
        #expand("workflow/out/filter/HouseholdTransmission-Stool-{sample}-filtered.1.fastq.gz",sample=samples),
        #expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/species/species_profile.txt",sample=samples),
        #"workflow/out/midasOutput/species/species_profile_all.txt"
        #expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),
        #expand("workflow/out/midasOutput/species/abundantSpecies_{household}.txt", household=households)
        #expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/snps/summary.txt", sample=samples)
        #expand("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_freq.txt.gz", species=snpsSpecies)
        #expand("workflow/report/calculateFixedDifferences/{species}/done.txt", species=snpsSpecies),
        #"workflow/report/calculateFixedDifferences/Bacteroides_stercoris_56735/done.txt"
        #expand("workflow/report/performStrainFishing/{species}/done.txt", species=snpsSpecies)
        #expand("workflow/report/calculateFixedDifferences/{species}/fixedDiffs_filterCoverage.txt.gz", species=snpsSpecies)
        #"workflow/report/performStrainFishing/Bacteroides_stercoris_56735/done.txt"
        #"workflow/report/performStrainFishing/Akkermansia_muciniphila_55290/done.txt"
        #"workflow/report/performStrainFishing/removeEmptyStrainFishing.done"
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-001.bam",
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-058.bam",
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-519.bam"
        #expand("workflow/out/T6SS/HouseholdTransmission-Stool-{sample}.bam", sample=samples),
        #"workflow/out/T6SS/T6SS-coverage.txt",
        #join(config["T6SSdir"],"T6SS-totalReads.txt")
        #join(config["MGEsdir"],"MGEs-totalReads.txt")
        #"workflow/out/MGEs/HouseholdTransmission-Stool-XBA-058.bam",
        #"workflow/out/MGEs/HouseholdTransmission-Stool-XBA-519.bam"
        #expand("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/XB.rds", species=XBspeciesOfInterest)
        #expand("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/lowerCoverage/{species}/snps_freq.txt.gz", species=snpsSpeciesLowerCoverage)
        #expand("workflow/report/performStrainFishing-lowerCoverage/{species}/done.txt", species=snpsSpeciesLowerCoverage)
        #expand("workflow/report/calculateFixedDifferences-lowerCoverage/{species}/done.txt", species=snpsSpeciesLowerCoverage)
        #"workflow/report/aggregateSNPsummaries/SNPsummaries.txt"
        #expand("workflow/report/performStrainFishing-lowerCoverage-coreGenesOnly/{species}/done.txt", species=snpsSpeciesLowerCoverage),
        #expand("workflow/report/calculateFixedDifferences-lowerCoverage-coreGenesOnly/{species}/done.txt", species=snpsSpeciesLowerCoverage)
        #"workflow/report/readDepths.txt"
        #expand("workflow/out/genomes/{species}/samples.txt", species=XBspeciesOfInterestGenomes),
        #expand("workflow/out/genomes/{species}/export.done.txt", species=XBspeciesOfInterestGenomes)
        #"workflow/out/genomes/genomeList.txt"
        #expand("{genomes}", genomes=XBgenomesamples)
        #expand("workflow/report/performStrainFishing-lowerCoverage-noSharedGenes/{species}/done.txt", species=snpsSpeciesLowerCoverage)
        expand("workflow/report/calculateFixedDifferences-lowerCoverage-noSharedGenes/{species}/done.txt", species=snpsSpeciesLowerCoverage)
        
        


#include: "workflow/rules/processRawReads.smk"
include: "workflow/rules/runMIDAS.smk"
include: "workflow/rules/postMIDAS.smk"
#include: "workflow/rules/mapReadsT6SS.smk"
#include: "workflow/rules/mapReadsMGEs.smk"
include: "workflow/rules/exportGenomesForAlignment.smk"
