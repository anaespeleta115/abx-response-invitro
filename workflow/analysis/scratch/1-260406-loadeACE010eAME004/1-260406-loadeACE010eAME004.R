library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library(patchwork)
library("phyloseq")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/out/"

# Read in full data file for eACE010 and eAME004
SR2601_data <- read.table("~/Documents/GitHub/abx-response-invitro/data/eACE010-eAME004ps_all.txt.gz", header = TRUE)



# Divide the data into corresponding experiments, eAME004, and add more metadata columns
eAME004_data <- SR2601_data %>% 
  filter(round2plate == "A", round1index %in% 2:7) %>% 
  mutate(replicate = case_when(
    str_detect(sample, "AME0060") ~ 1,
    str_detect(sample, "AME0061") ~ 2, 
    str_detect(sample, "AME0062") ~ 1,
    str_detect(sample, "AME0063") ~ 2,
    str_detect(sample, "AME0064") ~ 1,
    str_detect(sample, "AME0065") ~ 2
  )) %>% 
  rename(mixture = metadata) 


# Add recipient, donor, and day columns to eAME004 data
eAME004_data <- eAME004_data %>% 
  mutate(biosample2 = ifelse(str_detect(mixture, "blank") | str_detect(mixture, "B-mix"), 
                           str_sub(mixture, 1, 5),
                           str_sub(mixture, 1, 7)
         ),
  biosample1 = ifelse(str_detect(mixture, "blank") | str_detect(mixture, "B-mix"),
                          str_sub(mixture, 7, -1),
                          str_sub(mixture, 9, -1)
        ),
  day = str_sub(biosample1, -3),
  subject = str_sub(biosample1, 1, 3))
  

# eACE010
eACE010_data <- SR2601_data %>% 
  filter(round2plate %in% c("A", "B"), round1index %in% c(0,1)) %>% 
  rename(biosample = metadata) %>% 
  mutate(subject = str_sub(biosample, 1, 3),
         day = str_sub(biosample, -3))

# eACE010
eACE010_data_unpooled <- SR2601_data %>% 
  filter(round2plate %in% c("A"), round1index %in% c(0,1)) %>% 
  rename(biosample = metadata) %>% 
  mutate(subject = str_sub(biosample, 1, 3),
         day = str_sub(biosample, -3))


# eACE010
eACE010_data_pooled <- SR2601_data %>% 
  filter(round2plate %in% c("B"), round1index %in% c(0,1)) %>% 
  rename(biosample = metadata) %>% 
  mutate(subject = str_sub(biosample, 1, 3),
         day = str_sub(biosample, -3))


# adjust eACE010 days to be the correct timepoint
eACE010_data <- eACE010_data %>% 
  mutate(day = case_when(
    str_detect(day, "001") | str_detect(day, "002") | str_detect(day, "003") | str_detect(day, "022")| str_detect(day, "008") ~ "001",
    str_detect(day, "029") | str_detect(day, "028") | str_detect(day, "027") ~ "029",
    str_detect(day, "036") | str_detect(day, "037") | str_detect(day, "038") ~ "036",
    str_detect(day, "064")| str_detect(day, "063") | str_detect(day, "072") | str_detect(day, "059")| str_detect(day, "065") | str_detect(day, "066") ~ "064",
    TRUE ~ "0"
  ))


# Save full dataset into a text file

# eACE010 unpooled communities
write.table(eACE010_data_unpooled,
            file = paste0(OUTDIR,"eACE010_full_unpooled.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# eACE010 pooled communities
write.table(eACE010_data_pooled,
            file = paste0(OUTDIR,"eACE010_full_pooled.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)



# Define global variables
curr_replicate <- 1
limit_of_detection <-  1e-3
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 

# Filter the data to include only the selected replicate
eAME004_data <- eAME004_data %>% 
  filter(replicate == curr_replicate)
