
source("C:/abx-response-invitro/analysis/scratch/091125-prepareFAVAdata/091125-prepareFAVAdata.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/091125-runFAVAdonors/out/"

donor_subset <- combined %>% 
  mutate(biosample2 = ifelse(is.na(mixture), biosample1, biosample2)) %>% 
  filter(is.na(mixture) | str_sub(mixture, -3) == "029", biosample2 != "blank")

# Extract taxa-only matrix

taxa_cols <- intersect(colnames(donor_subset), taxa_names(phyloseq_obj))

taxa_mat <- as.matrix(donor_subset[, taxa_cols])


# Run FAVA on this subset

fava_val <- FAVA::fava(taxa_mat)
fava_val_norm <- FAVA::fava_norm(taxa_mat)

fava_val
fava_val_norm


# Run on all donors
fava_results <- donor_subset %>%
  group_by(biosample2) %>%
  summarise(
    fava      = FAVA::fava(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj)))))),
    fava_norm = FAVA::fava_norm(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj))))))
  )

# save into a text file
write_delim(fava_results, paste0(OUTDIR, "fava_donor_results_day29.txt"))

# save into a csv file
write.csv(
  fava_results,
  "C:/abx-response-invitro/analysis/scratch/091125-runFAVAdonors/out/fava_results_day29.csv",
  row.names = FALSE
)


# --------------- PLOT FAVA RESULTS -----------------------------------

fava_plot <- ggplot(fava_results, aes(x = 1, y = fava)) +
  geom_boxplot(fill = "steelblue") +
  # geom_jitter()+
  # geom_text(aes(label = biosample2), 
  #           position = position_jitter(height = 0), 
  #           vjust = -0.5, size = 1) +
  geom_text_repel(aes(label = biosample2),
    position = position_jitter(width = 0.1, height = 0),
    size = 2, segment.color = NA
  ) +
  labs(x = "", y = "FAVA") +
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x  = element_blank(),
                axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "fava_boxplot_day29"), fava_plot, 3, 3)


