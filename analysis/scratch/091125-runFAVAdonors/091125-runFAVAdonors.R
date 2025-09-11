
source("C:/abx-response-invitro/analysis/scratch/091125-prepareFAVAdata/091125-prepareFAVAdata.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/091125-runFAVAdonors/out/"


# # Select donor community of interest
# donor_subset <- combined %>%
#   filter(biosample1 == "XHB-029" | biosample2 == "XHB-029", is.na(mixture) | str_sub(mixture, -3) == "036", biosample2 != "blank")

donor_subset <- combined %>% 
  mutate(biosample2 = ifelse(is.na(mixture), biosample1, biosample2)) %>% 
  filter(is.na(mixture) | str_sub(mixture, -3) == "036", biosample2 != "blank")

# Extract taxa-only matrix

taxa_cols <- intersect(colnames(donor_subset), taxa_names(phyloseq_obj))

taxa_mat <- as.matrix(donor_subset[, taxa_cols])


# Run FAVA on this subset

fava_val <- FAVA::fava(taxa_mat)
fava_val_norm <- FAVA::fava_norm(taxa_mat)

fava_val
fava_val_norm


# run on all donors
fava_results <- donor_subset %>%
  group_by(biosample2) %>%
  summarise(
    fava      = FAVA::fava(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj)))))),
    fava_norm = FAVA::fava_norm(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj))))))
  )

# save into a text file
write_delim(fava_results, paste0(OUTDIR, "fava_donor_results.txt"))


