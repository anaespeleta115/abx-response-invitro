library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
library(rstatix)
library(PNWColors)

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

lastingResponses <- c("XBA", "XDA", "XEA", "XKA")

# Extract subject, day, household, and antibiotic information
e0026 <- e0026 %>%
  mutate(
    subject = str_sub(biosample1, 1, -5),
    day = str_sub(biosample1, -2),
    household = str_sub(biosample1, 1, -6),
    antibiotic = if_else(str_sub(biosample1, 1, -5) %in% subjectsAbx, 1, 0)
  )


# Divide up e0026 dataset into separate day datasets
e0026_day1 <- e0026 %>%  filter(str_detect(day, "01") | str_detect(day, "02") | str_detect(day, "03") | str_detect(day, "22")| str_detect(day, "08")) %>% 
  mutate(day = "01")
e0026_day29 <- e0026 %>%   filter(str_detect(day, "29") | str_detect(day, "28") | str_detect(day, "27")) %>% 
  mutate(day = "29")

e0026_day36 <- e0026 %>%  filter(str_detect(day, "36") | str_detect(day, "37")) %>% 
  mutate(day = "36")

e0026_day64 <- e0026 %>%  filter(str_detect(day, "64")| str_detect(day, "63") | str_detect(day, "72") | str_detect(day, "59")| str_detect(day, "65")) %>% 
  mutate(day = "64")

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
  select(biosample1, OTU, passage) %>%
  group_by(biosample1, passage) %>%
  dplyr::summarize(species_richness = n_distinct(OTU)) %>%
  left_join(
    combined_day_data %>% distinct(biosample1, day, subject, household, antibiotic),
    by = "biosample1"
  )


# Extract the top 25 families by relative abundance to make plots better
top_families <- e0026_day29 %>%
  group_by(Family) %>%
  dplyr::summarise(total_abundance = sum(relAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  # slice_head(n = 25) %>%
  pull(Family)

# Define color palette
my_colors <- readRDS("C:/abx-response-invitro/data/familyColorPalette.rds") 



