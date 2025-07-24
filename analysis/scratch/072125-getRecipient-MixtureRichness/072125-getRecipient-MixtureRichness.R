# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")




# Create recipient richness dataset
recipient_richness <- recipient_ASVs %>%
  select(biosample1, OTU) %>%
  distinct() %>%
  dplyr::count(biosample1, name = "species_richness")  %>% # count the distinct OTUs
  mutate(community = "recipient", 
         day = as.integer(str_extract(biosample1, "\\d{2}$")),
         subject = str_sub(biosample1, 1, -5)) %>%
  rename(sample = biosample1)


donor_richness <- single_donor_ASVs %>%
  select(biosample1, OTU) %>%
  distinct() %>%
  dplyr::count(biosample1, name = "species_richness")  %>% # count the distinct OTUs
  mutate(community = "donor", 
         day = as.integer(str_extract(biosample1, "\\d{2}$"))) %>%
  rename(sample = biosample1)


# Create mixture richness dataset
mixture_richness <- mixture_ASVs %>%
  filter(biosample2 %in% single_donor_ASVs$biosample1) %>%
  select(mixture, OTU) %>%
  distinct() %>%
  dplyr::count(mixture, name = "species_richness") %>% #count the distinct OTUs
  mutate(community = "mixture", 
         day = as.integer(str_extract(mixture, "\\d{2}$"))) %>%
  rename(sample = mixture)


# # Join richness datasets together to plot
# combined_richness <- bind_rows(recipient_richness, mixture_richness)







