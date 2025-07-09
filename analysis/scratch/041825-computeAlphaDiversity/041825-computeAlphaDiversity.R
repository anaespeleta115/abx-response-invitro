### Using Phyloseq to calculate diversity statistics


# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")


data <- readRDS("C:/abx-response-invitro/data/ps_all.rds")

hp_data_clean <- prune_taxa(taxa_sums(data) > 0, data)

e0026_obj <- subset_samples(hp_data_clean, experiment == "e0026")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/041825-computeAlphaDiversity/out/"

### Get alpha diversity metrics (ONLY RUN ONCE)

# Calculate alpha diversity.
calculateAlphaDiversity <- function(data) {

  # Use the estimate_richness function to calculate
  # an array of alpha-diversity statistics.
  alphaRaw <- estimate_richness(data, split=TRUE)
  # Add the sample names.
  alphaRaw$sample <- rownames(alphaRaw)
  # Calculate the Shannon effective number of species.
  alphaRaw <- alphaRaw %>%
	mutate(ShannonEffectiveSpecies=exp(Shannon))
  # Tidy the dataframe.
  alpha <- alphaRaw %>%
	mutate(sample=gsub("\\.","-",sample)) %>%
	pivot_longer(-sample, names_to="alphaStat", values_to="value")
}

alpha_diversity <- calculateAlphaDiversity(e0026_obj)
write_delim(alpha_diversity, paste0(OUTDIR, "alpha_diversity.txt"))




