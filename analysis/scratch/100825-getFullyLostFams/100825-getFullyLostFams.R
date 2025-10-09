# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")

curr_replicate <- 1

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


lost_family_abundance <- family_abundance_29 %>% 
  left_join(family_abundance_36, by = c("subject1", "Family")) %>% 
  mutate(fully_lost = ifelse(is.na(fam_abundance_36), 1, 0))

fully_lost_families <- lost_family_abundance %>% 
  filter(fully_lost == 1)





