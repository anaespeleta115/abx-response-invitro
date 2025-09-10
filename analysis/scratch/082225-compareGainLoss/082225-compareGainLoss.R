# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/082225-compareGainLoss/out/"

curr_replicate <- 1

## Calculate gain from mixture communities

# get relative abundance by family
fam_abundance_colonization <- actual_colonizers_results %>% 
  filter(day == "036"| day == "029", replicate == curr_replicate) %>% 
  group_by(day, Family, subject1, subject2) %>% 
  summarise(fam_relAbundance = sum(relAbundance))

## Check to make sure fam relative abundance adds to 1 in each mixture
# fam_abundance_colonization %>% 
#   group_by(day, subject1, subject2) %>% 
#   summarise(total_abundance = sum(fam_relAbundance))

# ----------------- Compute Post-Abx Family Gain ------------------------

gain_colonization <- fam_abundance_colonization %>% 
  pivot_wider(
    names_from = day,
    values_from = fam_relAbundance,
    names_prefix = "mix_fam_relAbundance_day_"
  )%>%
  replace_na(list(mix_fam_relAbundance_day_029 = 0, mix_fam_relAbundance_day_036 = 0)) %>%
  mutate(
    mix_fam_relAbundance_day_029 = if_else(mix_fam_relAbundance_day_029 == 0, 1e-4, mix_fam_relAbundance_day_029),
    mix_fam_relAbundance_day_036 = if_else(mix_fam_relAbundance_day_036 == 0, 1e-4, mix_fam_relAbundance_day_036),
    mix_ratio = mix_fam_relAbundance_day_036 / mix_fam_relAbundance_day_029,
    mix_log10_ratio = log10(mix_ratio)
  )


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_loss <- recipient_ASVs %>% 
  filter(day == "029" | day == "036", replicate == curr_replicate) %>% 
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  group_by(day, Family, subject1) %>% 
  summarise(fam_relAbundance = sum(relAbundance))

## Check for correct relAbundance calculation
# recipient_fam_loss %>%
#   group_by(biosample1) %>%
#   summarise(total_abundance = sum(fam_relAbundance))


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
    recipient_ratio = recipient_fam_relAbundance_day_036 / recipient_fam_relAbundance_day_029,
    recipient_log10_ratio = -log10(recipient_ratio)
  )


# ----------------- Combine gain and loss data ------------------------

compare_gain_loss <- gain_colonization %>% 
  left_join(recipient_fam_loss, by = c("Family", "subject1")) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(is.na(recipient_fam_relAbundance_day_029), 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(is.na(recipient_fam_relAbundance_day_036), 1e-4, recipient_fam_relAbundance_day_036),
    recipient_ratio = if_else(is.na(recipient_ratio), 1e-4, recipient_ratio),
    recipient_log10_ratio = if_else(is.na(recipient_log10_ratio), 1e-4, recipient_log10_ratio),
  )

# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

gain_loss_avg <- gain_colonization %>% 
  group_by(Family, subject1) %>% 
  summarise(avg_mix_ratio = mean(mix_log10_ratio)) %>% 
  left_join(recipient_fam_loss, by = c("Family", "subject1")) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(is.na(recipient_fam_relAbundance_day_029), 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(is.na(recipient_fam_relAbundance_day_036), 1e-4, recipient_fam_relAbundance_day_036),
    recipient_ratio = if_else(is.na(recipient_ratio), 1e-4, recipient_ratio),
    recipient_log10_ratio = if_else(is.na(recipient_log10_ratio), 1e-4, recipient_log10_ratio),
  ) %>% 
  filter(recipient_log10_ratio != 1e-04)

gain_loss_avg_plot <- gain_loss_avg %>% 
  # filter(recipient_fam_relAbundance_day_029 != 1e-04, recipient_fam_relAbundance_day_036 == 1e-04) %>% 
  ggplot(aes(x = recipient_log10_ratio, y = avg_mix_ratio, color = Family))+
  geom_point(size = 2)+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")


savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot"), gain_loss_avg_plot, 5, 6)



# ----------------- Plot x vs y ------------------------

# try filtering for families that were present in the recipient on day 29 and then completely lost on day 36

gain_loss_plot <- compare_gain_loss %>% 
  filter(recipient_fam_relAbundance_day_029 != 1e-04, recipient_fam_relAbundance_day_036 == 1e-04) %>%
  ggplot(aes(x = recipient_log10_ratio, y = mix_log10_ratio, color = Family))+
  geom_point(size = 2)+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")


savePNGPDF(paste0(OUTDIR, "gain_loss_plot_filtered"), gain_loss_plot, 5, 6)


  