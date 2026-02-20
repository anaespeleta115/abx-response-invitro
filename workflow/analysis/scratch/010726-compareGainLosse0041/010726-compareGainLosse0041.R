# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-gete0041colonization/102125-gete0041colonization.R")



# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/010726-compareGainLosse0041/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_abundance_pre_abx <- e0041_control_recipients %>%
  filter(recipient == "pre-abx", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(Family, recipient) %>%
  summarise(fam_relAbundance_pre_abx = sum(relAbundance))

# V1 family abundance
recipient_fam_abundance_V1 <- e0041_control_recipients %>%
  filter(recipient == "post-abx-V1", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(Family, recipient) %>%
  summarise(fam_relAbundance_V1 = sum(relAbundance))

#V2 family abundance
recipient_fam_abundance_V2 <- e0041_control_recipients %>%
  filter(recipient == "post-abx-V2", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(Family, recipient) %>%
  summarise(fam_relAbundance_V2 = sum(relAbundance))

#V1 family loss
recipient_fam_loss_V1 <- recipient_fam_abundance_pre_abx %>%
  left_join(recipient_fam_abundance_V1 %>% select(Family, fam_relAbundance_V1), by=c("Family")) %>% 
  replace_na(list(fam_relAbundance_V1 = 1e-4)) %>%
  mutate(
    recipient_ratio = fam_relAbundance_pre_abx / fam_relAbundance_V1 # either flip the fraction, or negate the ratio
  )

#V2 family loss
recipient_fam_loss_V2 <- recipient_fam_abundance_pre_abx %>%
  left_join(recipient_fam_abundance_V2 %>% select(Family, fam_relAbundance_V2), by=c("Family")) %>% 
  replace_na(list(fam_relAbundance_V2 = 1e-4)) %>%
  mutate(
    recipient_ratio = fam_relAbundance_pre_abx / fam_relAbundance_V2 # either flip the fraction, or negate the ratio
  )


# ----------------- Compute Post-Abx Family Gain (recipient 36 -> mixture 36) -------------------------------

# get V1 mixture family abundance
fam_abundance_mixture_V1 <- mixture_colonization_full %>% 
  filter(recipient == "post-abx-V1", replicate == curr_replicate) %>% 
  group_by(Family, recipient, donor) %>% 
  summarise(fam_relAbundance_mixture = sum(relAbundance))

# get V2 mixture family abundance
fam_abundance_mixture_V2 <- mixture_colonization_full %>% 
  filter(recipient == "post-abx-V2", replicate == curr_replicate) %>% 
  group_by(Family, recipient, donor) %>% 
  summarise(fam_relAbundance_mixture = sum(relAbundance))


# V1 recipient family abundance
recipient_fam_abundance_V1 <- e0041_control_recipients %>%
  filter(recipient == "post-abx-V1", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(Family, recipient) %>%
  summarise(fam_relAbundance_V1 = sum(relAbundance))

#V2 recipient family abundance
recipient_fam_abundance_V2 <- e0041_control_recipients %>%
  filter(recipient == "post-abx-V2", replicate == curr_replicate, relAbundance > 1e-4) %>%
  group_by(Family, recipient) %>%
  summarise(fam_relAbundance_V2 = sum(relAbundance))

#V1 family abundance gain in mixture
gain_mixture_V1 <- fam_abundance_mixture_V1 %>%
  left_join(recipient_fam_abundance_V1 %>% select(Family, fam_relAbundance_V1), by = c("Family")) %>% 
  replace_na(list(fam_relAbundance_mixture = 0, fam_relAbundance_V1 = 0)) %>%
  mutate(
    fam_relAbundance_V1 = if_else(fam_relAbundance_V1 == 0, 1e-4, fam_relAbundance_V1),
    fam_relAbundance_mixture = if_else(fam_relAbundance_mixture == 0, 1e-4, fam_relAbundance_mixture),
    mixture_ratio = fam_relAbundance_mixture/fam_relAbundance_V1
  ) 
# compute average mixture ratio across donors
gain_mixture_V1_avg <- gain_mixture_V1 %>% 
  group_by(Family, recipient) %>% 
  summarize(avg_mix_ratio = mean(mixture_ratio))

#V2 family abundance gain in mixture
gain_mixture_V2 <- fam_abundance_mixture_V2 %>%
  left_join(recipient_fam_abundance_V2 %>% select(Family, fam_relAbundance_V2), by = c("Family")) %>% 
  replace_na(list(fam_relAbundance_mixture = 0, fam_relAbundance_V2 = 0)) %>%
  mutate(
    fam_relAbundance_V2 = if_else(fam_relAbundance_V2 == 0, 1e-4, fam_relAbundance_V2),
    fam_relAbundance_mixture = if_else(fam_relAbundance_mixture == 0, 1e-4, fam_relAbundance_mixture),
    mixture_ratio = fam_relAbundance_mixture/fam_relAbundance_V2
  ) 
# compute average mixture ratio across donors
gain_mixture_V2_avg <- gain_mixture_V2 %>% 
  group_by(Family, recipient) %>% 
  summarize(avg_mix_ratio = mean(mixture_ratio))

# ----------------- Combine gain and loss data --------------------------------------------------------------
# Here, we left-join the gain data into the loss data because we only care about the families that were already present in the recipient
# V1 comparison
compare_gain_loss_V1 <- recipient_fam_loss_V1 %>%
  left_join(gain_mixture_V1 %>% select(-recipient, -fam_relAbundance_V1), by = c("Family")) %>%
  replace_na(list(fam_relAbundance_pre_abx = 1e-4)) %>%
  mutate(
    recipient_ratio = fam_relAbundance_pre_abx/fam_relAbundance_V1
  )

#V2 comparison
compare_gain_loss_V2 <- recipient_fam_loss_V2 %>%
  left_join(gain_mixture_V2 %>% select(-recipient, -fam_relAbundance_V2), by = c("Family")) %>%
  replace_na(list(fam_relAbundance_pre_abx = 1e-4)) %>%
  mutate(
    recipient_ratio = fam_relAbundance_pre_abx/fam_relAbundance_V2
  )


# # V1 comparison using full_join (W/ UNEXPECTED COLONIZERS)
# compare_gain_loss_V1 <- recipient_fam_loss_V1 %>%
#   full_join(gain_mixture_V1_avg %>% select(-recipient), by = c("Family")) %>%
#   replace_na(list(fam_relAbundance_pre_abx = 1e-4, recipient = "post_abx_V1", fam_relAbundance_V1 = 1e-4)) %>%
#   mutate(
#     recipient_ratio = fam_relAbundance_pre_abx/fam_relAbundance_V1
#   )
# 
# #V2 comparison using full_join (W/ UNEXPECTED COLONIZERS)
# compare_gain_loss_V2 <- recipient_fam_loss_V2 %>%
#   full_join(gain_mixture_V2_avg %>% select(-recipient), by = c("Family")) %>%
#   replace_na(list(fam_relAbundance_pre_abx = 1e-4, recipient = "post_abx_V2", fam_relAbundance_V2 = 1e-4)) %>%
#   mutate(
#     recipient_ratio = fam_relAbundance_pre_abx/fam_relAbundance_V2
#   )



# ----------------- Average gain ratio across mixtures --------------------------------------------------------------

# V1 PLOT -----------------
compare_gain_loss_V1 <- compare_gain_loss_V1 %>% 
  group_by(Family, recipient_ratio) %>%
  summarise(avg_mix_ratio = mean(mixture_ratio)) %>%
  filter(Family != "Hydrogenoanaerobacterium") #filter out this family because it was not present in any donor

gain_loss_avg_plot_V1 <- compare_gain_loss_V1 %>% 
  ggplot(aes(x = log10(recipient_ratio), y = log10(avg_mix_ratio), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  # labs(x = "log(pre-abx recipient fam relAbundance/recipient 36 fam relAbundance)", y = "log(mixture 36 fam relAbundance / recipient 36 fam relAbundance)")+
  labs(x = "Family loss in recipient 
       (log abundance)", y = "Average family gain in V1 mixture 
       (log abundance)", title = "Post-abx V1")+
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot_e0041_V1_attempt"), gain_loss_avg_plot_V1, 1.5, 2)



#V2 PLOT -----------------
compare_gain_loss_V2 <- compare_gain_loss_V2 %>%
  group_by(Family, recipient_ratio) %>%
  summarise(avg_mix_ratio = mean(mixture_ratio)) %>%
  filter(Family != "Hydrogenoanaerobacterium") #filter out this family because it was not present in any donor

gain_loss_avg_plot_V2 <- compare_gain_loss_V2 %>% 
  ggplot(aes(x = log10(recipient_ratio), y = log10(avg_mix_ratio), color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  # labs(x = "log(pre-abx recipient fam relAbundance/recipient 36 fam relAbundance)", y = "log(mixture 36 fam relAbundance / recipient 36 fam relAbundance)")+
  labs(x = "Family loss in recipient 
       (log abundance)", y = "Average family gain 
       in V2 mixture 
       (log abundance)", title = "Post-abx-V2")+
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot_e0041_V2_attempt"), gain_loss_avg_plot_V2, 1.5, 2)
