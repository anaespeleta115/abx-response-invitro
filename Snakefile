from os.path import join
import pandas as pd
import os

configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_table(config['sampleFileRaw'], header=None)
df.columns=["sample","samplelane","read1","read2","trim1","trim2"]

# Parse subject and household from sample names.
df['samplename']=df['sample']
df[['subject','timepoint']]=df.samplename.str.split("-",expand=True)
df['household']=df['subject'].astype(str).str[0:2]
df['studyArm']=df['subject'].astype(str).str[0:1]

# Generate lists of subjects and households.
subjects=list(set(df['subject'].tolist()))
households=list(set(df['household'].tolist()))
samplelanes=list(set(df['samplelane'].tolist()))
samples=list(set(df['samplename'].tolist()))

# Parse the list of species analyzed in the MIDAS snps module.
if(os.path.isdir("workflow/out/midasOutput/snps")):
    # Iterate through the species analyzed by the snps module.
    dirs=[x[0] for x in os.walk("workflow/out/midasOutput/snps/HouseholdTransmission-Stool")]
    # Remove the first element, which is the directory itself without subdirectories.
    dirs.pop(0)
    # Parse the species names and generate a list.
    snpsSpecies=[]
    for species in dirs:
        if "lowerCoverage" not in species:
            snpsSpecies.append(species.split("/")[5])
            
# Parse the list of species analyzed in the MIDAS snps module at lower coverage.
if(os.path.isdir("workflow/out/midasOutput/snps")):
    # Iterate through the species analyzed by the snps module.
    dirsLowerCoverage=[x[0] for x in os.walk("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/lowerCoverage")]
    # Remove the first element, which is the directory itself without subdirectories.
    dirsLowerCoverage.pop(0)
    # Parse the species names and generate a list.
    snpsSpeciesLowerCoverage=[]
    for species in dirsLowerCoverage:
        if "lowerCoverage" not in species.split("/")[6]:
            snpsSpeciesLowerCoverage.append(species.split("/")[6])
        
# Import the list of species of special interest in household XB.
dfXB=pd.read_table(config['XBspeciesOfInterest'], header=None)
dfXB.columns=["species"]
XBspeciesOfInterest=list(set(dfXB['species'].tolist()))

# Import the list of species whose genomes are of special interest in household XB.
dfXBgenomes=pd.read_table(config['XBspeciesOfInterestGenomes'], header=None)
dfXBgenomes.columns=["species"]
XBspeciesOfInterestGenomes=list(set(dfXBgenomes['species'].tolist()))

# Import the list of samples and species to export genomes for.
try:
    dfXBgenomesamples=pd.read_table("workflow/out/genomes/genomeList.txt", header=None)
    dfXBgenomesamples.columns=["genome"]
    XBgenomesamples=list(set(dfXBgenomesamples['genome'].tolist()))
except FileNotFoundError:
    print("workflow/out/genomes/genomeList.txtworkflow/out/genomes/genomeList.txt does not exist.")

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
