# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-gainLossAME004/out/"


# ----------------- Compute day 29 recipient family abundances (recipient 29) ------------------------

recipient_fam_abundance_29 <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(day, Family, subject) %>%
  summarise(fam_relAbundance = sum(relAbundance))

# ----------------- Compute Post-Abx family abundances (mixture 64) -------------------------------

recipient_fam_abundance_64 <- recipient_ASVs %>%
  filter(day == "064", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(day, Family, subject) %>%
  summarise(fam_relAbundance = sum(relAbundance))


# ----------------- Put both family abundance datasets together ------------------------

compare_recipient_mixture <- recipient_fam_abundance_29 %>% 
  left_join(recipient_fam_abundance_64, by = c("subject", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(fam_relAbundance.y = 1e-4, fam_relAbundance.x = 1e-4)) %>%
  # mutate(
  #   fam_relAbundance.y = if_else(fam_relAbundance.y == 0, 1e-4, fam_relAbundance.y),
  #   fam_relAbundance.x = if_else(fam_relAbundance.x == 0, 1e-4, fam_relAbundance.x)
  # ) %>%
  rename(recipient_abundance29 = fam_relAbundance.x, recipient_abundance64 = fam_relAbundance.y) %>% 
  filter(recipient_abundance29 != 1e-4) # filter out all taxa not present in day 29 recipient

# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

recipient_corr_avg_plot <- compare_recipient_mixture %>% 
  ggplot(aes(x = log(recipient_abundance29), y = log(recipient_abundance64), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Pre-abx recipient log family abundance", y = "One month post-abx log family abundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

recipient_corr_avg_plot

savePNGPDF(paste0(OUTDIR, "recipient_corr_avg"), recipient_corr_avg_plot, 5, 5)


