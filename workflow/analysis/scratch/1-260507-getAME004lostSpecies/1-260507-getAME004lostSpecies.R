source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260507-getAME004lostSpecies/out/"


# Make sure data is above limit of detection
eAME004_data <- eAME004_data %>% 
  filter(relAbundance > limit_of_detection)

recipients_day29 <- recipient_ASVs %>%
  filter(day == "029") %>%
  select(subject, OTU, day) %>%
  distinct()

recipients_day36 <- recipient_ASVs %>%
  filter(day == "036") %>%
  select(subject, OTU, day) %>%
  distinct() %>%
  mutate(present_day36 = 1)

# extract day 29 taxa not in that same recipient's day 36 taxa
lost_strains_36 <- recipients_day29 %>%
  left_join(recipients_day36, by = c("subject", "OTU")) %>%
  filter(is.na(present_day36)) %>%
  rename(day = day.x) %>% 
  select(subject, OTU, day) %>%
  mutate(lost_strain_36 = 1)

lost_strains_totals <- lost_strains_36 %>% 
  group_by(subject) %>% 
  summarize(total_species_lost = sum(lost_strain_36))

# Join lost strains back into recipient ASVs
recipient_ASVs <- recipient_ASVs %>% 
  left_join(lost_strains_36, by = c("subject", "day", "OTU")) %>%
  mutate(lost_strain_36 = ifelse(is.na(lost_strain_36), 0, lost_strain_36)) 

# Get lost taxa only
recipient_lost_36 <- recipient_ASVs %>%
  filter(lost_strain_36 == 1)


# Add lost tag to original dataset, plot compositions
AME004_lost_compositions <- eAME004_data %>% 
  filter(replicate == curr_replicate, biosample2 == "blank", day == "029") %>% 
  left_join(
    lost_strains_36 %>% select(subject, OTU, lost_strain_36),
    by = c("subject", "OTU")
  ) %>%
  mutate(lost_36 = ifelse(is.na(lost_strain_36), 0, 1)) %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = Family, alpha = factor(lost_36))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 1, "1" = 0)) +
  labs(x = "recipient", y = "relative abundance") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT


AME004_lost_compositions 

savePNGPDF(paste0(OUTDIR, "eAME004_maintained_compositions"), AME004_lost_compositions, 3, 4)


# Plot compositions of differential colonizers
AME004_diff_compositions <- colonization_results %>% 
  filter(day == "036") %>% 
  ggplot(aes(x = biosample2, y = relAbundance, fill = Family, alpha = factor(colonized_day36)))+
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(x = "day 36 mixture", y = "relative abundance")+
  facet_wrap(~biosample1)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

AME004_diff_compositions
  


# put both recipient and donor compositions together
eAME004_joint_compositions <- AME004_lost_compositions / AME004_diff_compositions + plot_layout(heights = c(1, 3))

savePNGPDF(paste0(OUTDIR, "eAME004_lost-colonizers_compositions"), eAME004_joint_compositions, 6, 4)

