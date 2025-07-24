# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")



recipient_ASVs <- recipient_ASVs %>%
  mutate(day = str_extract(biosample1, "\\d{3}"))

recipients_day29 <- recipient_ASVs %>%
  filter(day == "029", replicate == 1) %>%
  select(biosample1, OTU) %>%
  distinct()

recipients_day36 <- recipient_ASVs %>%
  filter(day == "036", replicate == 1) %>%
  select(biosample1, OTU) %>%
  distinct()

recipients_day64 <- recipient_ASVs %>%
  filter(day == "064", replicate == 1) %>%
  select(biosample1, OTU) %>%
  distinct()

# extract day 29 taxa not in day 36 taxa
lost_strains_29_36 <- recipients_day29 %>% 
  filter(!(OTU %in% recipients_day36$OTU)) %>% 
  mutate(lost_strain_29_36 = 1)

# extract day 29 taxa not in day 64 taxa
lost_strains_29_64 <- recipients_day29 %>% 
  filter(!(OTU %in% recipients_day64$OTU)) %>% 
  mutate(lost_strain_29_64 = 1)

# extract day 36 taxa not in day 64 taxa
lost_strains_36_64 <- recipients_day36 %>% 
  filter(!(OTU %in% recipients_day64$OTU)) %>% 
  mutate(lost_strain_36_64 = 1)

# Join in lost strains from each of the day categories back into recipient ASVs
recipient_ASVs <- recipient_ASVs %>%
  left_join(lost_strains_29_36, by = c("biosample1", "OTU")) %>%
  mutate(lost_strain_29_36 = ifelse(is.na(lost_strain_29_36), 0, lost_strain_29_36))

recipient_ASVs <- recipient_ASVs %>%
  left_join(lost_strains_29_64, by = c("biosample1", "OTU")) %>%
  mutate(lost_strain_29_64 = ifelse(is.na(lost_strain_29_64), 0, lost_strain_29_64))

recipient_ASVs <- recipient_ASVs %>%
  left_join(lost_strains_36_64, by = c("biosample1", "OTU")) %>%
  mutate(lost_strain_36_64 = ifelse(is.na(lost_strain_36_64), 0, lost_strain_36_64))


# Get lost taxa only
recipient_lost_29_36 <- recipient_ASVs %>%
  filter(lost_strain_29_36 == 1)

recipient_lost_29_64 <- recipient_ASVs %>%
  filter(lost_strain_29_64 == 1)

recipient_lost_36_64 <- recipient_ASVs %>%
  filter(lost_strain_36_64 == 1)
