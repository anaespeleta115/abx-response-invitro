library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
library(PNWColors)
library(paletteer)


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/010326-loade0026MetagenomicData/out"


metaGdata <- read.table("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loadMetaGenomes/out/e0026_metaG_full.txt", sep = "\t", header = TRUE, quote = "", comment.char = "")

# Specify antibiotic-taking subjects
subjectsAbx <- c("XAA","XBA","XCA","XDA","XEA","XFA","XGA","XHB","XIA","XJA",
                 "XKA","XLA","XMA","XNA","XOA","XPA","XQA","XRA","XSA","XTA",
                 "XUA","XVA")



# Set limit of detection
limit_of_detection <- 1e-3

# Add subject and antibiotic columns
metaGdata <- metaGdata %>% 
  mutate(subject = str_sub(sample, 1, 3),
         antibiotic = ifelse(subject %in% subjectsAbx, 1, 0),
         day = str_sub(sample, -3))

# Adjust days to be the correct timepoint
metaGdata <- metaGdata %>% 
  mutate(day = case_when(
    str_detect(day, "001") | str_detect(day, "002") | str_detect(day, "003") | str_detect(day, "022")| str_detect(day, "008") ~ "001",
    str_detect(day, "029") | str_detect(day, "028") | str_detect(day, "027") ~ "029",
    str_detect(day, "036") | str_detect(day, "037") | str_detect(day, "038") ~ "036",
    str_detect(day, "064")| str_detect(day, "063") | str_detect(day, "072") | str_detect(day, "059")| str_detect(day, "065") | str_detect(day, "066") ~ "064",
    TRUE ~ "0"
  ))

# Set metagenomic color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/commonFamiliesPalette.rds") 

# filter out low abundance taxa
metaGdata_filtered <- metaGdata %>% 
  filter(relative_abundance > limit_of_detection)
