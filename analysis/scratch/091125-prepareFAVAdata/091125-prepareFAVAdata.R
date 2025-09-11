# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")

library(phyloseq)
library(FAVA)


# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/091125-prepareFAVAdata/out/"


num_replicate <- 1


# Get phyloseq objects
phyloseq_obj <- readRDS("C:/abx-response-invitro/data/ps_all.rds")


# Donor-only samples (from e0026)
ps_donors <- subset_samples(phyloseq_obj, experiment == "e0026" & passage == 8 & biosample1 %in% donor_communities)

# Mixtures (from e0029)
ps_mixtures <- subset_samples(phyloseq_obj, experiment == "e0029" & replicate == num_replicate)


donor_relab <- FAVA::relab_phyloseq(ps_donors)
mixt_relab  <- FAVA::relab_phyloseq(ps_mixtures)

# Match donors to mixtures
# Separate metadata and taxa
donor_taxa <- donor_relab %>%
  mutate(type = "donor") %>%
  select(biosample1, everything())

mixt_taxa <- mixt_relab %>%
  mutate(type = "mixture") %>%
  select(mixture, everything())

# Join donors to their mixtures
combined <- bind_rows(donor_relab, mixt_relab)



