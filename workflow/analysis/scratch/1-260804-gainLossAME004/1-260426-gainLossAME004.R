# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260804-gainLossAME004/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_loss <- recipient_ASVs %>%
  filter(day == "029" | day == "036", replicate == curr_replicate, relAbundance > 1e-4) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(day, Family, subject) %>%
  summarise(fam_relAbundance = sum(relAbundance))

recipient_fam_loss <- recipient_fam_loss %>%
  pivot_wider(
    names_from = day,
    values_from = fam_relAbundance,
    names_prefix = "recipient_fam_relAbundance_day_"
  )%>%
  replace_na(list(recipient_fam_relAbundance_day_029 = 0, recipient_fam_relAbundance_day_036 = 0)) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(recipient_fam_relAbundance_day_029 == 0, 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(recipient_fam_relAbundance_day_036 == 0, 1e-4, recipient_fam_relAbundance_day_036),
    recipient_ratio = recipient_fam_relAbundance_day_029 / recipient_fam_relAbundance_day_036,
    recipient_log10_ratio = log10(recipient_ratio)
  )

# ----------------- Compute Post-Abx Family Gain (recipient 36 -> mixture 36) -------------------------------

fam_abundance_mixture <- colonization_results %>% 
  mutate(subject2 = str_sub(biosample2, 1, 3)) %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  group_by(day, Family, subject, subject2) %>% 
  summarise(fam_relAbundance = sum(relAbundance)) 

recipient_fam_36 <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>%
  group_by(day, Family, subject) %>%
  summarise(fam_relAbundance = sum(relAbundance)) 

gain_mixture <- recipient_fam_36 %>%
  full_join(fam_abundance_mixture, by = c("subject", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(fam_relAbundance.y = 1e-4, fam_relAbundance.x = 1e-4)) %>%
  rename(recipient_abundance36 = fam_relAbundance.x, mixture_abundance36 = fam_relAbundance.y) %>% 
  mutate(
    mix_ratio = mixture_abundance36/recipient_abundance36, # y = mix 36, x = recipient 36
    mix_log10_ratio = log10(mix_ratio)
  ) 


# ----------------- Combine gain and loss data --------------------------------------------------------------

compare_gain_loss <- gain_mixture %>% 
  full_join(recipient_fam_loss, by = c("Family", "subject")) %>% 
  mutate(
    recipient_ratio = if_else(is.na(recipient_ratio), 1e-4, recipient_ratio),
    recipient_log10_ratio = if_else(is.na(recipient_log10_ratio), 1e-4, recipient_log10_ratio),
  ) %>% 
  filter()


# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

gain_loss_avg <- gain_mixture %>% 
  group_by(Family, subject) %>% 
  summarise(avg_mix_ratio = mean(mix_log10_ratio)) %>% 
  left_join(recipient_fam_loss, by = c("Family", "subject")) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(is.na(recipient_fam_relAbundance_day_029), 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(is.na(recipient_fam_relAbundance_day_036), 1e-4, recipient_fam_relAbundance_day_036),
    recipient_ratio = if_else(is.na(recipient_ratio), 1e-4, recipient_ratio),
    recipient_log10_ratio = if_else(is.na(recipient_log10_ratio), 1e-4, recipient_log10_ratio),
  ) 

gain_loss_avg_plot <- gain_loss_avg %>%
  ggplot(aes(x = recipient_log10_ratio, y = avg_mix_ratio, color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Recipient loss (log(recipient day 29 relAbundance / 
       recipient day 36 relAbundance)", y = "Average mixture gain (log(mixture day 36 relAbundance / 
       recipient 36 relAbundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

gain_loss_avg_plot

savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot"), gain_loss_avg_plot, 3, 4)



# ----------------- Plot x vs y -------------------------------------------------------------------------------

# try filtering for families that were present in the recipient on day 29 and then completely lost on day 36

gain_loss_plot <- compare_gain_loss %>% 
  # filter(recipient_fam_relAbundance_day_029 != 1e-04, recipient_fam_relAbundance_day_036 == 1e-04) %>%
  ggplot(aes(x = recipient_log10_ratio, y = mix_log10_ratio, color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  labs(x = "Recipient loss (log(recipient day 29 relAbundance / 
       recipient day 36 relAbundance)", y = "Mixture gain (log(mixture day 36 relAbundance / 
       recipient 36 relAbundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT



savePNGPDF(paste0(OUTDIR, "gain_loss_plot"), gain_loss_plot, 3, 4)


# ----------------- Classify loss  -------------------------------------------------------------------------------


compare_gain_loss <- compare_gain_loss %>%
  mutate(
    category = case_when(
      recipient_log10_ratio > 0.5 ~ "Lost",
      recipient_log10_ratio < -0.5 ~ "Gained",
      TRUE ~ "Stable"
    )
  )

# Summarize average mixture gain per category
compare_gain_loss %>%
  group_by(subject1, subject2, category) %>%
  summarise(
    mean_gain = mean(mix_log10_ratio, na.rm = TRUE),
    median_gain = median(mix_log10_ratio, na.rm = TRUE),
    n = n()
  )

