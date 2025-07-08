library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)



household_data <- read.table("C:/Users/anaes/OneDrive/UCI_Spring25/rotation/e0026-e0029-e0030.txt", header = TRUE)


# Divide dataset into separate tables by experiment. Filter out unnecessary columns
e0026 <- filter(household_data, str_detect(household_data$sample, "e0026")) %>% 
  select(sample, biosample1, experiment, passage, OTU, count, replicate, relAbundance, Phylum, Family, Genus)


# List the antibiotic-taking subjects.
# Based on metagenomic and metabolomic data, annotate XHB as antibiotic-taking
# and XHC as a control subject.
subjectsAbx <- c("XAA","XBA","XCA","XDA","XEA","XFA","XGA","XHB","XIA","XJA",
                 "XKA","XLA","XMA","XNA","XOA","XPA","XQA","XRA","XSA","XTA",
                 "XUA","XVA")

# Extract subject, day, household, and antibiotic information
e0026 <- e0026 %>%
  mutate(
    subject = str_sub(biosample1, 1, -5),
    day = str_sub(biosample1, -3),
    household = str_sub(biosample1, 1, -6),
    antibiotic = if_else(str_sub(biosample1, 1, -5) %in% subjectsAbx, 1, 0)
  )


# Divide up e0026 dataset into separate day datasets
e0026_day1 <- e0026 %>%  filter(str_detect(biosample1, "001") | str_detect(biosample1, "002") | str_detect(biosample1, "003") | str_detect(biosample1, "022")| str_detect(biosample1, "008")) %>% 
  mutate(day = "001")

e0026_day29 <- e0026 %>%   filter(str_detect(biosample1, "029") | str_detect(biosample1, "028") | str_detect(biosample1, "027")) %>% 
  mutate(day = "029")

e0026_day36 <- e0026 %>%  filter(str_detect(biosample1, "036") | str_detect(biosample1, "037")) %>% 
  mutate(day = "036")

e0026_day64 <- e0026 %>%  filter(str_detect(biosample1, "064")| str_detect(biosample1, "063") | str_detect(biosample1, "072") | str_detect(biosample1, "059")| str_detect(biosample1, "065")) %>% 
  mutate(day = "064")

# Combine all day datasets
combined_day_data <- bind_rows(e0026_day1, e0026_day29, e0026_day36, e0026_day64)

# Get all samples that have every passage sequenced
e0026_all_passages <- e0026 %>%
  group_by(biosample1)  %>% 
  mutate(num_passages = n_distinct(passage)) %>% 
  filter(num_passages == 9) %>% 
  ungroup() 

# Compute species richness
e0026_richness <- combined_day_data %>%
  filter(relAbundance > 0.001) %>%
  select(biosample1, OTU) %>%
  distinct() %>%
  dplyr::count(biosample1, name = "species_richness") %>%
  left_join(
    combined_day_data %>% distinct(biosample1, day, subject, household, antibiotic),
    by = "biosample1"
  )

# Extract the top 25 families by relative abundance to make plots better
top_families <- e0026_day29 %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(relAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 25) %>%
  pull(Family)

# Define color palette
my_colors <- readRDS("C:/abx-response-invitro/data/familyColorPalette.rds") 



