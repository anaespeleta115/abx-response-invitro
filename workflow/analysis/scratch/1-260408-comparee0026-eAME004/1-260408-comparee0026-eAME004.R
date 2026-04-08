source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/040325-loadData/loadData.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-comparee0026-eAME004/out/"


# Keep only passage 8 and replicate 1 of e0026.
e0026_p8_r1 <- e0026 %>% 
  filter(passage == 8, replicate == 1) %>% 
  mutate(experiment == "e0026") %>%  # add experiment column
  select(sample, biosample1, experiment, OTU, count, relAbundance, Phylum, Family, Genus, subject, day, -antibiotic)

# Keep only relevant columns from eACE010
eACE010_subset <- eACE010_data %>% 
  mutate(experiment = "eACE010") %>% 
  rename(biosample1 = biosample) %>% 
  select(sample, biosample1, experiment, OTU, count, relAbundance, Phylum, Family, Genus, subject, day)

# Bind the two datasets together
e0026_eACE010_data <- rbind(eACE010_subset, e0026_p8_r1)



# Group by sample and family
summarized_e0026_eACE010_data <- e0026_eACE010_data %>% 
  filter(day == "029") %>% 
  group_by(experiment, subject, day, Family) %>% 
  summarize(total_abundance = sum(relAbundance), .groups = "drop")



# Summarize relAbundance per OTU per subject per day
fc_e0026_eACE010 <- summarized_e0026_eACE010_data %>%
  pivot_wider(
    names_from = experiment,
    values_from = total_abundance,
    names_prefix = "relAbundance_exp_"
  ) %>%
  replace_na(list(relAbundance_exp_e0026 = 1e-4, relAbundance_exp_eACE010 = 1e-4)) %>% 
  mutate(fold_change = relAbundance_exp_eACE010/relAbundance_exp_e0026) %>% 
  mutate(log10_fc = log10(fold_change))


# Plot fold change
p_fc_e0026_eACE010 <- fc_e0026_eACE010 %>%
  ggplot(aes(x = fct_reorder(Family, -log10_fc), y = log10_fc, fill = Family)) +
  geom_boxplot(width = 0.9, outlier.shape = NA) +
  scale_fill_manual(values = my_colors)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred")+ 
  labs(
    title = "",
    x = "Family",
    y = "Family Log10 Fold Change e0026/eACE010"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "e0026_eACE010_family_fold_change"), p_fc_e0026_eACE010, 4, 6)
