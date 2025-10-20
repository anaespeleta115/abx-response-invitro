# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/092525-compareLossGainDifference/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

# recipient_fam_loss <- recipient_ASVs %>%
#   filter(day == "029" | day == "036", replicate == curr_replicate, relAbundance > 1e-4) %>%
#   mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
#   group_by(day, Family, subject1) %>%
#   summarise(fam_relAbundance = sum(relAbundance))

recipient_fam_abundance_29 <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate, relAbundance > 1e-4) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(Family, subject1) %>%
  summarise(fam_relAbundance_29 = sum(relAbundance))

recipient_fam_abundance_36 <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate, relAbundance > 1e-4) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(Family, subject1) %>%
  summarise(fam_relAbundance_36 = sum(relAbundance))


recipient_fam_loss <- recipient_fam_abundance_29 %>%
  left_join(recipient_fam_abundance_36, by=c("subject1", "Family")) %>% 
  replace_na(list(fam_relAbundance_36 = 1e-4)) %>%
  mutate(
    log_recipient_abundance_29 = log10(fam_relAbundance_29),
    log_recipient_abundance_36 = log10(fam_relAbundance_36),
    recipient_ratio = fam_relAbundance_29 / fam_relAbundance_36, # either flip the fraction, or negate the ratio
  )

# recipient_fam_loss <- recipient_fam_loss %>%
#   pivot_wider(
#     names_from = day,
#     values_from = fam_relAbundance,
#     names_prefix = "recipient_fam_relAbundance_day_"
#   )%>%
#   replace_na(list(recipient_fam_relAbundance_day_029 = 0, recipient_fam_relAbundance_day_036 = 0)) %>%
#   mutate(
#     recipient_fam_relAbundance_day_029 = if_else(recipient_fam_relAbundance_day_029 == 0, 1e-4, recipient_fam_relAbundance_day_029),
#     recipient_fam_relAbundance_day_036 = if_else(recipient_fam_relAbundance_day_036 == 0, 1e-4, recipient_fam_relAbundance_day_036),
#     recipient_ratio = recipient_fam_relAbundance_day_029 / recipient_fam_relAbundance_day_036   # check this !!!
#     # recipient_diff = log10(recipient_fam_relAbundance_day_029)
#   )


# ----------------- Compute Post-Abx Family Gain (recipient 36 -> mixture 36) -------------------------------

fam_abundance_mixture <- actual_colonizers_results %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  group_by(day, Family, subject1, subject2) %>% 
  summarise(fam_relAbundance_mixture = sum(relAbundance)) %>% 
  rename(day_mixture = day)

recipient_fam_36 <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(day, Family, subject1) %>%
  summarise(fam_relAbundance_recipient_36 = sum(relAbundance))  %>% 
  rename(day_recipient = day)

gain_mixture <- fam_abundance_mixture %>%
  left_join(recipient_fam_36, by = c("subject1", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(fam_relAbundance_mixture = 0, fam_relAbundance_recipient_36 = 0, day_recipient = "036")) %>%
  mutate(
    fam_relAbundance_recipient_36 = if_else(fam_relAbundance_recipient_36 == 0, 1e-4, fam_relAbundance_recipient_36),
    log_mix_abundance_36 = log10(fam_relAbundance_mixture),
    log_recipient_abundance_36 = log10(fam_relAbundance_recipient_36),
    mixture_ratio = fam_relAbundance_mixture/fam_relAbundance_recipient_36
  ) 


# ----------------- Combine gain and loss data --------------------------------------------------------------

compare_gain_loss <- gain_mixture %>% 
  left_join(recipient_fam_loss %>%  dplyr::select(-fam_relAbundance_36), by = c("Family", "subject1", "log_recipient_abundance_36")) %>% 
  replace_na(list(fam_relAbundance_29 = 1e-4)) %>% 
  mutate(
    log_recipient_abundance_29 = log10(fam_relAbundance_29),
    recipient_ratio = fam_relAbundance_29/fam_relAbundance_recipient_36
  )
  # filter(!(log_recipient_abundance_29 == -4.0 & log_recipient_abundance_36 == -4.0)) # filter out the families that weren't in the recipient neither before 
# after the antibiotic and that were gained in the mixture., optionally also filter out families not present on day 29 in the recipient

# ----------------- (SIDE QUEST: get unexpected colonizers) --------------------------------------------------------------

unexpected_colonizers <- gain_mixture %>% 
  left_join(recipient_fam_loss %>%  dplyr::select(-fam_relAbundance_36), by = c("Family", "subject1", "log_recipient_abundance_36")) %>% 
  replace_na(list(fam_relAbundance_29 = 1e-4)) %>% 
  mutate(
    log_recipient_abundance_29 = log10(fam_relAbundance_29),
    recipient_ratio = fam_relAbundance_29/fam_relAbundance_recipient_36
  ) %>% 
  filter(log_recipient_abundance_29 == -4.0 & log_recipient_abundance_36 == -4.0)


unexpected_colonizers_plot <- unexpected_colonizers %>% 
  ggplot(aes(x = fct_reorder(Family, - fam_relAbundance_mixture), y = fam_relAbundance_mixture))+
  geom_col(fill = "pink", color = "black")+
  labs(x = "Family",
       y = "Day 36 mixture relAbundance of family")+
  facet_grid(~subject1)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))

savePNGPDF(paste0(OUTDIR, "unexpected_colonizers"), unexpected_colonizers_plot, 3, 4)

  
# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

gain_loss_avg <- compare_gain_loss %>% 
  group_by(Family, subject1, recipient_ratio) %>% 
  summarise(avg_mix_ratio = mean(mixture_ratio))
  
gain_loss_avg_plot <- gain_loss_avg %>% 
  ggplot(aes(x = log(recipient_ratio), y = log(avg_mix_ratio), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  # labs(x = "log(recipient 29 fam relAbundance/recipient 36 fam relAbundance)", y = "log(mixture 36 fam relAbundance / recipient 36 fam relAbundance)")+
  labs(x = "Family loss in day 36 recipient 
       (log abundance)", y = "Average family gain in day 36 mixture 
       (log abundance)")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot"), gain_loss_avg_plot, 3, 4)


