# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")

library(phyloseq)
library(FAVA)


# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/090925-FAVAanalysis/out/"


num_replicate <- 1


# Get phyloseq objects
phyloseq_obj <- readRDS("C:/abx-response-invitro/data/ps_all.rds")


e0026_obj <- subset_samples(
  phyloseq_obj,
  passage == 8 & experiment == "e0026" & replicate == num_replicate
)

e0026_otu_table <- otu_table(e0026_obj)

e0029_obj <- subset_samples(
  phyloseq_obj,
  experiment == "e0029" & replicate == num_replicate
)

# Transform abundances
e0026_data <- FAVA::relab_phyloseq(e0026_obj) %>% 
  select(experiment, sample, mixture, 20:ncol(e0026_data)) 


e0029_data <- FAVA::relab_phyloseq(e0029_obj) %>% 
  select(experiment, sample, mixture, 20:ncol(e0029_data)) 

# 1) Subset your ps to the cohort you want
e0029_obj <- subset_samples(phyloseq_obj, experiment == "e0029" & replicate == num_replicate & biosample1 == "XEA-036" & biosample2 != "blank")


# 2) Convert to relative abundance table with metadata on the left
relab_e0029 <- FAVA::relab_phyloseq(e0029_obj)  # data.frame: metadata cols + taxa cols
# Taxa columns = those matching phyloseq taxa_names
taxa_cols_e0029 <- intersect(colnames(relab_e0029), taxa_names(phyloseq_obj))
# Metadata + taxa
metadata_cols_e0029 <- setdiff(colnames(relab_e0029), taxa_cols_e0029)
# 3) Select only the taxa columns (drop metadata)
relab_metadata_e0029 <- relab_e0029 %>% select(all_of(metadata_cols_e0029))
relab_taxa_e0029 <- relab_e0029 %>% select(all_of(taxa_cols_e0029))
taxa_mat_e0029 <- as.matrix(relab_taxa_e0029)  # rows = samples, cols = taxa; rows sum to 1



# 4) Run FAVA
fava_val <- FAVA::fava(taxa_mat_e0029)
fava_val_norm <- FAVA::fava_norm(taxa_mat_e0029)





# # adding rows to the object:
# sample_data(phyloseq_obj)$NewMetric <- external_vector
# names(external_vector) <- sample_names(phyloseq_obj)




