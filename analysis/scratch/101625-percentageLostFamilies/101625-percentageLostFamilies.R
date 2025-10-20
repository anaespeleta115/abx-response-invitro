# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/101625-percentageLostFamilies/out/"

curr_replicate <- 1

# ------------------------------------------ % abundance lost out of total family abundance -----------------------------------------------

recipient_ASVs <- recipient_ASVs %>%
  filter(replicate == curr_replicate) %>% 
  mutate(day = str_extract(biosample1, "\\d{3}"))

recipients_day29_families <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate) %>% #  add replicate 1 filter
  select(biosample1, Family, OTU, relAbundance) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  distinct()

recipients_day36_families <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>% #  add replicate 1 filter
  select(biosample1, Family, OTU, relAbundance) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  distinct()

family_abundance_29 <- recipients_day29_families %>% 
  group_by(subject1, Family) %>% 
  summarize(fam_abundance_29 = sum(relAbundance))


family_abundance_36 <- recipients_day36_families %>% 
  group_by(subject1, Family) %>% 
  summarize(fam_abundance_36 = sum(relAbundance))


family_abundance <- family_abundance_29 %>% 
  left_join(family_abundance_36, by = c("subject1", "Family")) %>% 
  mutate(fam_abundance_36 = ifelse(is.na(fam_abundance_36), 0, fam_abundance_36), 
         difference = as.numeric(fam_abundance_29) - as.numeric(fam_abundance_36)) %>% 
  mutate(percentage_lost = difference/fam_abundance_29,
         lost = ifelse(percentage_lost > 0, 1, 0),
         gained = ifelse(percentage_lost < 0, 1, 0))

# ------------------------------------------ % of family lost (in abundance) -----------------------------------------------

lost_family_abundance <- family_abundance %>% 
  filter(lost == 1) 

percent_lost_family_plot <- lost_family_abundance %>% filter(subject1 == "XBA") %>% 
  select(subject1, Family, fam_abundance_29, percentage_lost) %>% 
  ggplot(aes(x = fct_reorder(Family, -fam_abundance_29), y = percentage_lost, fill = fam_abundance_29)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of day 29 family abundance lost",
       fill = "Day 29 family abundance")+
  facet_wrap(~subject1, scales = "free_y")+
  scale_fill_viridis_c(option = "magma")+ # limits = c(0.0001, 0.3)
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_lost_fams_day36_XBA"), percent_lost_family_plot , 4, 4)

# ------------------------------------------ % of family gained (in abundance) -----------------------------------------------

gained_family_abundance <- family_abundance %>% 
  filter(gained == 1)

# Use full join to capture those taxa not present on day 29, that were gained on day 36
family_abundance_full <- family_abundance_29 %>% 
  full_join(family_abundance_36, by = c("subject1", "Family")) %>% 
  mutate(fam_abundance_36 = ifelse(is.na(fam_abundance_36), 0, fam_abundance_36),
         fam_abundance_29 = ifelse(is.na(fam_abundance_29), 1e-4, fam_abundance_29),
         difference = as.numeric(fam_abundance_29) - as.numeric(fam_abundance_36)) %>% 
  mutate(percentage_lost = difference/fam_abundance_29,
         lost = ifelse(percentage_lost > 0, 1, 0),
         gained = ifelse(percentage_lost < 0, 1, 0))

gained_family_abundance <- family_abundance_full %>% 
  filter(gained == 1)

percentage_gained_family_plot <- gained_family_abundance %>% filter(subject1 == "XBA") %>% 
  select(subject1, Family, fam_abundance_29, percentage_lost) %>% 
  ggplot(aes(x = fct_reorder(Family, -fam_abundance_29), y = -percentage_lost, fill = fam_abundance_29)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of day 29 family abundance gained",
       fill = "Day 29 family abundance")+
  facet_wrap(~subject1, scales = "free_y")+
  scale_fill_viridis_c(option = "magma")+ # limits = c(0.0001, 0.3)
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_gained_fams_day36_XBA"), percentage_gained_family_plot , 4, 4)


# ------------------------------------------ % OTUs lost per family -----------------------------------------------

recipients_day29_OTU_counts <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate) %>%
  select(biosample1, Family, OTU) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  group_by(subject1, Family) %>% 
  summarize(OTUs_per_fam_29 = n())


recipients_day36_OTU_counts <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>%
  select(biosample1, Family, OTU) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  group_by(subject1, Family) %>% 
  summarize(OTUs_per_fam_36 = n())


family_OTU_counts <- recipients_day29_OTU_counts %>% 
  left_join(recipients_day36_OTU_counts, by = c("subject1", "Family")) %>% 
  mutate(OTUs_per_fam_36 = ifelse(is.na(OTUs_per_fam_36), 0, OTUs_per_fam_36), 
         difference = as.numeric(OTUs_per_fam_29) - as.numeric(OTUs_per_fam_36)) %>% 
  mutate(percentage_OTU_lost = difference/OTUs_per_fam_29,
         lost = ifelse(percentage_OTU_lost > 0, 1, 0),
         gained = ifelse(percentage_OTU_lost < 0, 1, 0))



hist_OTU_loss <- family_OTU_counts %>% filter(lost == 1) %>% 
  mutate(percentage_OTU_lost = percentage_OTU_lost * 100) %>%  # convert to %
  ggplot(aes(x = percentage_OTU_lost)) +
  geom_histogram(bins = 5, fill = "skyblue", color = "black", alpha = 0.7) +
  # geom_density(color = "darkblue", size = 1) +
  facet_wrap(~subject1) +
  labs(
    title = "",
    x = "% of OTUs Lost (Day 29 to Day 36)",
    y = "Number of families"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "OTU-loss-percentages-histogram"), hist_OTU_loss , 4, 4)

density_OTU_loss <- family_OTU_counts %>% filter(lost == 1) %>% 
  mutate(percentage_OTU_lost = percentage_OTU_lost * 100) %>%  # convert to %
  ggplot(aes(x = percentage_OTU_lost)) +
  # geom_histogram(bins = 10, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(~subject1) +
  labs(
    title = "",
    x = "% of OTUs Lost (Day 29 to Day 36)",
    y = "Density (relative frequency of families)"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "OTU-loss-percentages-density"), density_OTU_loss , 4, 4)

# Notes on this probability density: 

# A higher density at a given x means many families have loss values near that %.
# A lower density means fewer families have loss values near that %.
# The width of the curve captures the spread of those values across your subjects.
  
