# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/082725-compareDeltaRichnessColonization/out/"

curr_replicate <- 1



# ----------------- Mixture Colonization Gain --------------------------------

## get day 29 mixture family colonizer OTU counts
fam_colonizer_counts_mixture_29 <- actual_colonizers_results %>% 
  filter(day == "029", replicate == curr_replicate, actual_colonizer == 1) %>% 
  group_by(Family, subject1, subject2) %>% 
  summarise(fam_OTUs_29_mix = n_distinct(OTU))

## get day 36 mixture family colonizer OTU counts
fam_colonizer_counts_mixture_36 <- actual_colonizers_results %>% 
  filter(day == "036", replicate == curr_replicate, actual_colonizer == 1) %>% 
  group_by(Family, subject1, subject2) %>% 
  summarise(fam_OTUs_36_mix = n_distinct(OTU))

## join data from both days 
family_colonizer_gain_mixture <- fam_colonizer_counts_mixture_36 %>% 
  full_join(fam_colonizer_counts_mixture_29, by = c("Family", "subject1", "subject2")) %>%
  replace_na(list(fam_OTUs_29_mix = 0, fam_OTUs_36_mix = 0)) %>%
  mutate(
    delta_colonization_mixture = fam_OTUs_36_mix - fam_OTUs_29_mix
  )

# ------------------ Recipient Loss -----------------------------------------

## get day 29 recipient family OTU counts
fam_OTU_counts_recipient_29 <- recipient_ASVs %>% 
  filter(day == "029", replicate == curr_replicate) %>% 
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  group_by(Family, subject1) %>% 
  summarise(fam_OTUs_29_recipient = n_distinct(OTU))



## get day 36 recipient family OTU counts
fam_OTU_counts_recipient_36 <- recipient_ASVs %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  group_by(Family, subject1) %>% 
  summarise(fam_OTUs_36_recipient = n_distinct(OTU))


## join data from both days
family_OTU_loss_recipient <- fam_OTU_counts_recipient_36 %>% 
  full_join(fam_OTU_counts_recipient_29, by = c("Family", "subject1")) %>%
  replace_na(list(fam_OTUs_29_recipient = 0, fam_OTUs_36_recipient = 0)) %>%
  mutate(
    delta_OTU_recipient = fam_OTUs_29_recipient - fam_OTUs_36_recipient
  )

# -------------- compare Delta Loss and Colonization -------------------------------


deltaLossColonization <- family_colonizer_gain_mixture %>% 
  left_join(family_OTU_loss_recipient, by = c("Family", "subject1")) %>% 
  mutate(
    fam_OTUs_29_recipient = if_else(is.na(fam_OTUs_29_recipient), 0, fam_OTUs_29_recipient),
    fam_OTUs_36_recipient = if_else(is.na(fam_OTUs_36_recipient), 0, fam_OTUs_36_recipient),
    delta_OTU_recipient = if_else(is.na(delta_OTU_recipient), 0, delta_OTU_recipient)
  )


# Plot

delta_species_colonization_plot <- deltaLossColonization %>% 
  ggplot(aes(x = delta_OTU_recipient, y = delta_colonization_mixture, color = Family))+
  geom_point(size = 2)+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")


savePNGPDF(paste0(OUTDIR, "delta_species_colonization_plot"), delta_species_colonization_plot, 5, 6)

# --------------- compare delta loss vs. median colonization -------------------------

deltaLossColonization_avg <- deltaLossColonization %>% 
  group_by(Family, subject1) %>% 
  summarise(avg_colonization = mean(delta_colonization_mixture)) %>% 
  left_join(family_OTU_loss_recipient, by = c("Family", "subject1")) %>% 
  mutate(
    fam_OTUs_29_recipient = if_else(is.na(fam_OTUs_29_recipient), 0, fam_OTUs_29_recipient),
    fam_OTUs_36_recipient = if_else(is.na(fam_OTUs_36_recipient), 0, fam_OTUs_36_recipient),
    delta_OTU_recipient = if_else(is.na(delta_OTU_recipient), 0, delta_OTU_recipient)
  )


deltaLossColonization_avg_plot <- deltaLossColonization_avg %>% 
  # filter(recipient_fam_relAbundance_day_029 != 1e-04, recipient_fam_relAbundance_day_036 == 1e-04) %>% 
  ggplot(aes(x = delta_OTU_recipient, y = avg_colonization, color = Family))+
  geom_point(size = 1)+
  scale_color_manual(values = my_colors)+
  labs(x = "")+
  facet_wrap(~subject1)+
  theme(legend.position = "none")


savePNGPDF(paste0(OUTDIR, "deltaLossColonization_avg_plot"), deltaLossColonization_avg_plot, 3.5, 4)

  
