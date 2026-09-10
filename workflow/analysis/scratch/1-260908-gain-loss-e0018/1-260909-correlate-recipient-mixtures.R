source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260908-gain-loss-e0018/out/"



# ----------------- Compute Pre-abx Recipient Abundance (recipient 29) ------------------------

recipient_fam_abundance_preabx <- preabx_recipient_OTUs %>%
  group_by(Family, household) %>%
  summarise(recipient_fam_abundance = sum(relAbundance))

# ----------------- Compute Post-Abx Family Abundance (mixture 36) -------------------------------

fam_abundance_mixture <- mixture_OTUs %>%
  group_by(Family, household) %>%
  summarise(
    mixture_fam_abundance = sum(relAbundance),
    .groups = "drop"
  )

# ----------------- Put both family abundance datasets together ------------------------

compare_recipient_mixture <- recipient_fam_abundance_preabx %>% 
  left_join(fam_abundance_mixture, by = c("household", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(recipient_fam_abundance = 1e-4, mixture_fam_abundance = 1e-4))


# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

recipient_mixture_plot <- compare_recipient_mixture %>% 
  filter(household %in% c("L", "P", "A", "Q")) |> 
  ggplot(aes(x = log(recipient_fam_abundance), y = log(mixture_fam_abundance), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  # geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "log(Pre-abx recipient family abundance)", y = "log(Post-abx mixture family abundance)")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~household)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

recipient_mixture_plot

savePNGPDF(paste0(OUTDIR, "subset_log_recipient_mixture_plot"), recipient_mixture_plot, 3, 3)


