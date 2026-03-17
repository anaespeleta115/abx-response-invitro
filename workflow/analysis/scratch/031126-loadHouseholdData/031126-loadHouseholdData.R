library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
library(PNWColors)
library(paletteer)


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loadHouseholdData/out/"


household_data <- read.table("~/Documents/GitHub/abx-response-invitro/data/household_speciesAbundances.txt", header = TRUE)

# Specify antibiotic-taking subjects
subjectsAbx <- c("XAA","XBA","XCA","XDA","XEA","XFA","XGA","XHB","XIA","XJA",
                 "XKA","XLA","XMA","XNA","XOA","XPA","XQA","XRA","XSA","XTA",
                 "XUA","XVA")



# Set limit of detection
limit_of_detection <- 1e-3

# Add subject and antibiotic columns
household_data <- household_data %>% 
  mutate(subject = str_sub(sample, 1, 3),
         antibiotic = ifelse(subject %in% subjectsAbx, 1, 0),
         day = str_sub(sample, -3))


# Adjust all the timepoints that are off
household_data <- household_data %>%
  mutate(day = case_when(
    str_detect(day, "001") | str_detect(day, "002") | str_detect(day, "003") ~ "001",
    str_detect(day, "029") ~ "029",
    str_detect(sample, "XTB-028") ~ "029",
    str_detect(sample, "XTB-038") ~ "036",
    str_detect(day, "036") ~ "036",
    str_detect(day, "064")| str_detect(day, "063") | str_detect(day, "072") | str_detect(day, "059")| str_detect(day, "065") ~ "064",
    TRUE ~ "0"
  )) %>% 
  filter(sample != "XRA-059")

# Set metagenomic color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/commonFamiliesPalette.rds") 

# filter out low abundance taxa
household_filtered <- household_data %>% 
  filter(relative_abundance > limit_of_detection, str_detect(sample, "X"))


