source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-geteAME004Richness/1-260408-geteAME004SpeciesRichness.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-plotRichnessColonization/out/"



# adjust recipient diversity df
eAME004_recipient_diversity <- eAME004_recipient_diversity %>% 
  mutate(day = str_sub(day, -2, -1), day = paste0("0", day))

# get total colonization
total_colonization <- colonization_results %>%
  group_by(mixture, biosample1) %>%
  summarize(total_colonizers = sum(colonizer)) %>% 
  mutate(
    subject = str_sub(biosample1, 1, -5),
    subject2 = str_sub(mixture, 1, 3),
    day = str_sub(biosample1, -3))


# average colonization for each subject
avg_colonization <- total_colonization %>%
  group_by(biosample1, subject, day) %>% 
  summarize(avg_colonization = mean(total_colonizers))

# join together recipient diversity and colonization
diversity_colonization_df <- avg_colonization %>%
  left_join(eAME004_recipient_diversity, by = c("subject", "day"))
  

# ------------ PLOT ----------------


diversity_colonization_plot <- diversity_colonization_df %>% 
  ggplot(aes(x = nASVs, y = avg_colonization, color = factor(day))) +
  geom_point(size = 2)+
  geom_text_repel(aes(label = biosample1), size = 2, max.overlaps = 20)+
  labs(
    x = "Recipient species richness",
    y = "Avg number of colonizers",
    color = "Study day"
  )+
  # facet_wrap(~subject)+
  scale_color_discrete()+
  DEFAULTS.THEME_PRINT +
  theme(legend.key.height = unit(7, "pt"),
        legend.key.spacing = unit(2, "pt"))


diversity_colonization_plot

savePNGPDF(paste0(OUTDIR, "eAME004_colonization-richness_avg_no-facets"), diversity_colonization_plot, 4, 4)
