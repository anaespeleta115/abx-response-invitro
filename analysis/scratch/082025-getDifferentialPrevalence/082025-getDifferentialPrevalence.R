# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/082025-getDifferentialPrevalence/out/"

curr_replicate <- 1

# ------------------------------ DONOR ANALYSIS ----------------------------------------------

# Get family prevalence across donor communities
prevalence_donors <- single_donor_ASVs %>% 
  group_by(Family) %>%
  summarise(num_donors_present = n_distinct(biosample1))

# Determine prevalence threshold
donor_prevalence_plot <- prevalence_donors %>%
  ggplot(aes(x = fct_reorder(Family, -num_donors_present), y = num_donors_present)) +
  geom_col(fill = "lightblue", color = "black")+
  labs(x = "Family",
       y = "Number of donor communities family is present")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 7, angle = 90))

savePNGPDF(paste0(OUTDIR, "donor_prevalence_plot"), donor_prevalence_plot, 4, 6)


quantile(prevalence_donors$num_donors_present)
  
# Filter by that threshold (adjust later if necessary)
prevalence_donors <- prevalence_donors %>%
  filter(num_donors_present >= 6)

# Extract list of prevalent Families
prevalence_donor_list <- prevalence_donors$Family


abundance_donors <- single_donor_ASVs %>% 
  group_by(Family, biosample1) %>% 
  summarize(total_relAbundance = sum(relAbundance))

# Determine prevalence threshold
donor_abundance_plot <- abundance_donors %>% filter(biosample1 %in% c("XKB-029", "XJB-029", "XHB-029", "XIB-029")) %>% 
  ggplot(aes(x = fct_reorder(Family, -total_relAbundance), y = total_relAbundance, fill = Family)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "Total family abundance")+
  facet_wrap(~biosample1)+
  scale_fill_manual(values = my_colors)+
  DEFAULTS.THEME_PRINT+
  theme(    axis.text.x = element_blank(),   # removes text labels
            axis.ticks.x = element_blank(), legend.position = "none")

savePNGPDF(paste0(OUTDIR, "donor_abundance_plot_row3"), donor_abundance_plot, 3, 4)



# ------------------------------ DAY 36 ----------------------------------------------

# Get differential prevalence across mixtures (DAY 36)
differential_prevalence_mixes_36 <- actual_colonizers_results %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  group_by(Family, subject1, subject2) %>% 
  summarise(num_otus_differential = sum(diff_colonizer_36)) %>% 
  mutate(otus_differential = ifelse(num_otus_differential > 0, 1, 0)) %>% 
  summarise(num_mixes_differential = sum(otus_differential)) %>% 
  filter(Family %in% prevalence_donor_list)

total_diff_prevalence_36 <- prevalence_donors %>% 
  left_join(differential_prevalence_mixes_36, by = "Family") %>% 
  mutate(family_differential_prevalence = num_mixes_differential/num_donors_present)

# Get overall prevalence scores (out of the most prevalent families in the donor communities, 
# this is their mean differential prevalence across all mixtures)
mean_prevalence_36 <- total_diff_prevalence_36 %>% 
  group_by(Family) %>% 
  summarise(mean_prevalence = mean(family_differential_prevalence))


# Find family prevalence in the lost taxa
prevalence_recipients <- recipient_ASVs %>% 
  filter(day == "029", replicate == curr_replicate) %>% 
  group_by(Family) %>%
  summarise(num_recipients_present = n_distinct(biosample1))

# Determine prevalence threshold
recipient_prevalence_plot <- prevalence_recipients %>%
  ggplot(aes(x = fct_reorder(Family, -num_recipients_present), y = num_recipients_present)) +
  geom_col(fill = "lightpink", color = "black")+
  labs(x = "Family",
       y = "Number of recipient communities family is present")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 7, angle = 90))

savePNGPDF(paste0(OUTDIR, "recipient_prevalence_plot_29"), recipient_prevalence_plot, 4, 6)


# Filter by that threshold (adjust later if necessary)
prevalence_recipients <- prevalence_recipients %>%
  filter(num_recipients_present >= 3)

# Extract list of prevalent Families
prevalence_recipient_list <- prevalence_recipients$Family


# Get lost family prevalence across mixtures (DAY 36)
lost_families_36 <- recipient_ASVs %>% 
  filter(day == "029", replicate == curr_replicate) %>% 
  group_by(Family, biosample1) %>% 
  summarise(
    total_OTUs = n_distinct(OTU),
    num_otus_lost = sum(lost_strain_29_36),
    total_relAbundance = sum(relAbundance),
    lost_relAbundance = sum(relAbundance[lost_strain_29_36 == 1])
  ) %>% 
  mutate(prop_lost = num_otus_lost / total_OTUs, prop_relAbundance_lost = lost_relAbundance/total_relAbundance)

# Plot number of OTUs lost per family per recipient
lost_families_plot_36 <- lost_families_36 %>%
  filter(Family %in% prevalence_donor_list) %>% 
  ggplot(aes(x = fct_reorder(Family, - prop_lost), y = prop_lost)) +
  geom_col(fill = "lightyellow", color = "black")+
  labs(x = "Family",
       y = "Fraction of OTUs lost")+
  facet_wrap(~biosample1)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))

savePNGPDF(paste0(OUTDIR, "lost_families_plot_36"), lost_families_plot_36, 4, 5)

lost_family_relAbundance_prop_plot_36 <- lost_families_36 %>% 
  filter(Family %in% prevalence_donor_list) %>% 
  ggplot(aes(x = fct_reorder(Family, - total_relAbundance), y = prop_relAbundance_lost)) +
  geom_col(fill = "lightyellow", color = "black")+
  labs(x = "Family",
       y = "Fraction of family relAbundance lost")+
  facet_wrap(~biosample1)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))

savePNGPDF(paste0(OUTDIR, "lost_relAbundance_prop_plot_36"), lost_family_relAbundance_prop_plot_36 , 4, 5)


lost_family_relAbundance_plot_36 <- lost_families_36 %>% 
  filter(Family %in% prevalence_donor_list) %>% 
  ggplot(aes(x = fct_reorder(Family, - lost_relAbundance), y = lost_relAbundance)) +
  geom_col(fill = "thistle", color = "black")+
  labs(x = "Family",
       y = "RelAbundance of OTUs lost")+
  facet_wrap(~biosample1)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))

savePNGPDF(paste0(OUTDIR, "lost_relAbundance_plot_36"), lost_family_relAbundance_plot_36 , 4, 5)

# Get lost family prevalence
lost_prevalence_36 <- lost_families_36 %>%
  mutate(otus_lost = ifelse(num_otus_lost > 0, 1, 0)) %>%
  group_by(Family, biosample1) %>%
  summarise(num_communities_lost = sum(otus_lost)) %>%
  filter(Family %in% prevalence_recipient_list)


# Out of the total recipient communities the family was present in, how many 
total_lost_prevalence_36 <- prevalence_recipients %>%
  left_join(lost_prevalence_36, by = "Family") %>%
  mutate(family_lost_prevalence = num_communities_lost/num_recipients_present)


# Plot differential prevalence

diff_prevalence_plot_36 <- total_diff_prevalence_36 %>%
  ggplot(aes(x = fct_reorder(Family, -family_differential_prevalence), y = family_differential_prevalence)) +
  geom_col(fill = "lightgreen", color = "black")+
  labs(x = "Family",
       y = "Differential prevalence across recipient mixtures")+
  facet_wrap(~subject1)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 7, angle = 90))

savePNGPDF(paste0(OUTDIR, "diff_prevalence_plot_36"), diff_prevalence_plot_36, 4, 5)



# Plot lost prevalence


lost_prevalence_plot_36 <- total_lost_prevalence_36 %>%
  filter(Family %in% prevalence_donor_list) %>% 
  ggplot(aes(x = fct_reorder(Family, -family_lost_prevalence), y = family_lost_prevalence)) +
  geom_col(fill = "lightyellow", color = "black")+
  labs(x = "Family",
       y = "Lost family prevalence across recipient mixtures")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 7, angle = 90))

savePNGPDF(paste0(OUTDIR, "lost_prevalence_plot_36"), lost_prevalence_plot_36, 4, 5)



# # Get overall prevalence scores (out of the most prevalent families in the donor communities, 
# # this is their mean differential prevalence across all mixtures)
# mean_prevalence_36 <- total_prevalence_36 %>% 
#   group_by(Family) %>% 
#   summarise(mean_prevalence = mean(family_differential_prevalence))


# ------------------------------ DAY 64 ----------------------------------------------

# Get differential prevalence across mixtures (DAY 64)
differential_prevalence_mixes_64 <- actual_colonizers_results %>% 
  filter(day == "064", replicate == curr_replicate) %>% 
  group_by(Family, subject1, subject2) %>% 
  summarise(num_otus_differential = sum(diff_colonizer_64)) %>% 
  mutate(otus_differential = ifelse(num_otus_differential > 0, 1, 0)) %>% 
  summarise(num_mixes_differential = sum(otus_differential)) %>% 
  filter(Family %in% prevalence_donor_list)

total_prevalence_64 <- prevalence_mixes_64 %>% 
  left_join(differential_prevalence_mixes_64, by = "Family") %>% 
  mutate(family_differential_prevalence = num_mixes_differential/num_donors_present)


# Find family prevalence in the lost taxa


