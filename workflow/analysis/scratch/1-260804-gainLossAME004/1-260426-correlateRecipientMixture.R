# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260804-gainLossAME004/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_abundance_29 <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(day, Family, subject) %>%
  summarise(fam_relAbundance = sum(relAbundance))

# ----------------- Compute Post-Abx Family Gain (recipient 36 -> mixture 36) -------------------------------

fam_abundance_mixture <- colonization_results %>%
  mutate(subject2 = str_sub(biosample2, 1, 3)) %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  group_by(day, Family, subject, subject2) %>% 
  summarise(fam_relAbundance = sum(relAbundance)) 

# ----------------- Put both family abundance datasets together ------------------------

compare_recipient_mixture <- recipient_fam_abundance_29 %>% 
  left_join(fam_abundance_mixture, by = c("subject", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(fam_relAbundance.y = 1e-4, fam_relAbundance.x = 1e-4)) %>%
  rename(recipient_abundance29 = fam_relAbundance.x, mixture_abundance36 = fam_relAbundance.y)


avg_fam_abundance_mixture <- compare_recipient_mixture %>% 
  group_by(Family, subject, recipient_abundance29) %>% 
  summarize(avg_fam_abundance_mixture = mean(mixture_abundance36)) %>% 
  filter(recipient_abundance29 != 1e-4) # filter out all taxa not present in day 29 recipient


# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

recipient_mixture_avg_plot <- avg_fam_abundance_mixture %>% 
  ggplot(aes(x = log(recipient_abundance29), y = log(avg_fam_abundance_mixture), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Pre-abx recipient log family abundance", y = "Post-abx mixture log family abundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

recipient_mixture_avg_plot

savePNGPDF(paste0(OUTDIR, "recipient29-mixture36_avg_plot"), recipient_mixture_avg_plot, 5, 5)


