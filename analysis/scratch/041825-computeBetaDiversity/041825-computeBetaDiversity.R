# Get beta diversity metrics (ONLY RUN ONCE)

# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")


data <- readRDS("C:/abx-response-invitro/data/ps_all.rds")

hp_data_clean <- prune_taxa(taxa_sums(data) > 0, data)

e0026_obj <- subset_samples(hp_data_clean, experiment == "e0026")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/041825-computeBetaDiversity/out/"


# List the distance-calculation methods to be used.
distMethods <- c("jsd","bray","jaccard")

# Write a function to calculate a distance matrix using the specified method
# and convert the data into tidy format.
calculateBeta <- function(data, distMethod) {
  # Calculate the distance matrix using the specified method.
  betaRaw <- distance(data, method=distMethod)

  # Convert distance matrix to a dataframe.
  beta <- as.matrix(betaRaw)
  beta <- as.data.frame(beta)
  beta$sample1 <- rownames(beta)
  # Tidy the dataframe.
  beta <- beta %>%
	pivot_longer(-sample1, names_to="sample2", values_to="value")
  beta <- beta %>%
	filter(sample1 != sample2) %>%
	mutate(method=distMethod)
}

# Calculate the distance matrix for all of the specified methods on the species abundances.
# Combine the distance matrices for all methods.
betaSpecies <- do.call(rbind, lapply(distMethods, function(distMethod) {
  print(distMethod)
  calculateBeta(e0026_obj, distMethod)
}))
# Export the distance matrix generated for all of the sample pairs
# using all of the specified methods.
write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta.txt.gz"))



