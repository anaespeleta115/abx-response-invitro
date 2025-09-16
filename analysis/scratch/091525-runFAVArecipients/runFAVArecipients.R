
source("C:/abx-response-invitro/analysis/scratch/091125-prepareFAVAdata/091125-prepareFAVAdata.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/091525-runFAVArecipients/out/"



recipient_subset <- combined_recipients %>% 
  filter(biosample1 != "blank", str_sub(mixture, -3) == "029")

# Extract taxa-only matrix

taxa_cols <- intersect(colnames(recipient_subset), taxa_names(phyloseq_obj))

taxa_mat <- as.matrix(recipient_subset[, taxa_cols])



# Run on all recipients
fava_results <- recipient_subset %>%
  group_by(biosample1) %>%
  summarise(
    fava      = FAVA::fava(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj)))))),
    fava_norm = FAVA::fava_norm(as.matrix(across(all_of(intersect(colnames(.), taxa_names(phyloseq_obj))))))
  )


# save into a csv file
write.csv(
  fava_results,
  "C:/abx-response-invitro/analysis/scratch/091525-runFAVArecipients/out/fava_results_day36.csv",
  row.names = FALSE
)


# --------------- PLOT FAVA RESULTS -----------------------------------

fava_plot <- ggplot(fava_results, aes(x = 1, y = fava)) +
  geom_boxplot(fill = "#f8766d") +
  # geom_jitter()+
  # geom_text(aes(label = biosample2), 
  #           position = position_jitter(height = 0), 
  #           vjust = -0.5, size = 1) +
  geom_text_repel(aes(label = biosample1),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 3, segment.color = NA
  ) +
  labs(x = "", y = "FAVA") +
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "fava_boxplot_day29"), fava_plot, 3, 3)


