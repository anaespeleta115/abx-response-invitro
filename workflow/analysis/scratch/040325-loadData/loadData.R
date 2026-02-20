library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
# library(rstatix)
library(PNWColors)
library(paletteer)

# scale_colour_paletteer_d("lisa::BridgetRiley")
# scale_color_paletteer_d("lisa::BridgetRiley")
# scale_fill_paletteer_d("lisa::BridgetRiley")
# paletteer_d("lisa::BridgetRiley")

household_data <- read.table("~/Documents/GitHub/abx-response-invitro/data/e0026-e0029-e0030.txt", header = TRUE)

# List the antibiotic-taking subjects.
# Based on metagenomic and metabolomic data, annotate XHB as antibiotic-taking
# and XHC as a control subject.
subjectsAbx <- c("XAA","XBA","XCA","XDA","XEA","XFA","XGA","XHB","XIA","XJA",
                 "XKA","XLA","XMA","XNA","XOA","XPA","XQA","XRA","XSA","XTA",
                 "XUA","XVA")

# Divide dataset into separate tables by experiment. Filter out unnecessary columns
e0026 <- filter(household_data, str_detect(household_data$sample, "e0026")) %>% 
  select(sample, biosample1, experiment, passage, OTU, count, replicate, relAbundance, Phylum, Family, Genus) %>%
  mutate(subject = str_sub(biosample1, 1, 3),
         antibiotic = ifelse(subject %in% subjectsAbx, 1, 0))


write.table(
  e0026,
  file = "~/Documents/GitHub/abx-response-invitro/analysis/scratch/040325-loadData/e0026_subset.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


lastingResponses <- c("XBA", "XDA", "XEA", "XKA")

limit_of_detection <-  1e-3

# Extract subject, day, household, and antibiotic information
e0026 <- e0026 %>%
  filter(biosample1 != "blank") %>% 
  mutate(
    subject = str_sub(biosample1, 1, 3),
    day = str_sub(biosample1, -3),
    household = str_sub(biosample1, 1, -6),
    antibiotic = if_else(str_sub(biosample1, 1, -5) %in% subjectsAbx, 1, 0)
  )


e0026 <- e0026 %>% 
  mutate(day = case_when(
    str_detect(day, "001") | str_detect(day, "002") | str_detect(day, "003") | str_detect(day, "022")| str_detect(day, "008") ~ "001",
    str_detect(day, "029") | str_detect(day, "028") | str_detect(day, "027") ~ "029",
    str_detect(day, "036") | str_detect(day, "037") ~ "036",
    str_detect(day, "064")| str_detect(day, "063") | str_detect(day, "072") | str_detect(day, "059")| str_detect(day, "065") ~ "064",
    TRUE ~ "0"
  ))


# # Combine all day datasets
# combined_day_data <- bind_rows(e0026_day1, e0026_day29, e0026_day36, e0026_day64)

# Get all samples that have every passage sequenced
e0026_all_passages <- e0026 %>%
  group_by(biosample1)  %>% 
  mutate(num_passages = n_distinct(passage)) %>% 
  filter(num_passages == 9) %>% 
  ungroup() 


e0026_species_richness <- e0026 %>% 
  filter(relAbundance > limit_of_detection) %>% 
  group_by(passage, subject, day, antibiotic) %>% 
  summarize(species_richness = n())

# # Extract the top 25 families by relative abundance to make plots better
# top_families <- e0026_day29 %>%
#   group_by(Family) %>%
#   dplyr::summarise(total_abundance = sum(relAbundance, na.rm = TRUE)) %>%
#   arrange(desc(total_abundance)) %>%
#   # slice_head(n = 25) %>%
#   pull(Family)

# Define color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 



