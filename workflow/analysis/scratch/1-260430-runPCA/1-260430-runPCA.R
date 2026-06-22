source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260430-runPCA/out/"




# get updated metadata ready for this object
new_meta <- eAME004_data %>%
  select(sample, replicate, mixture) %>%
  distinct()

# extract old metadata from obj
meta <- data.frame(sample_data(eAME004_obj))
meta$sample <- rownames(meta)

# join them together by sample
meta_updated <- meta %>%
  left_join(new_meta, by = "sample") %>% 
  filter(round2plate == "A") %>% 
  mutate(subject = str_sub(metadata, -7, -5), 
         day = str_sub(metadata, -3, -1))


# keep only samples in meta_updated
eAME004_obj_filtered <- prune_samples(meta_updated$sample, eAME004_obj)

# make sample names the rownames again
rownames(meta_updated) <- meta_updated$sample
meta_updated$sample <- NULL

# update sample metadata in the filtered object
sample_data(eAME004_obj_filtered) <- sample_data(meta_updated)


## ------------------- Run PCA on phyloseq obj ---------------------


# distance + PCoA
ordu <- ordinate(eAME004_obj_filtered, method = "PCoA", distance = "jsd")


# plot with sample data
pca <- plot_ordination(eAME004_obj, ordu) +
  geom_point(size = 2) +
  geom_text(aes(label = metadata), vjust = -0.5, size = 3)+
  facet_wrap(~round1index)

savePNGPDF(paste0(OUTDIR, "eAME004-eACE010"), pca, 10, 10)


### --------------- Subset by round1index (only plot ACE data) --------------------

# subset_obj <- subset_samples(eAME004_obj_filtered, round1index %in% c("0", "1") & 
#                                day %in% c("028", "029", "036", "064") &
#                                metadata %in% c(donor_communities, recipient_ASVs$biosample1))


subset_obj <- subset_samples(eAME004_obj_filtered, round1index %in% c("0", "1", "2", "4", "6") &
                             metadata %in% c(donor_communities, recipient_ASVs$mixture))


subset_obj <- prune_taxa(taxa_sums(subset_obj) > 0, subset_obj)

ordu <- ordinate(subset_obj, method = "PCoA", distance = "jsd")

ACE_pca <- plot_ordination(subset_obj, ordu, color = "subject") +
  geom_point(size = 2) +
  geom_text_repel(aes(label = metadata), vjust = -0.5, size = 2)

savePNGPDF(paste0(OUTDIR, "eAME004_pca"), ACE_pca, 6, 6)



