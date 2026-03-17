library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
library(PNWColors)
library(paletteer)



# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loadMetaGenomes/out/"

# Get e0026 metagenomic data
metaGdata <- read.table("~/Documents/GitHub/abx-response-invitro/workflow/out/midasOutput/species/species_profile_all.txt", header = TRUE)

# Get midas taxonomy table
midas_taxonomy <- read.table("~/Downloads/species_taxonomy.txt", sep = "\t", header = TRUE, quote = "", comment.char = "")

midas_taxonomy <- midas_taxonomy %>% 
  mutate(family = ifelse(species_id == "Alistipes_onderdonkii_55464", "Rikenellaceae", family))

# select only relevant columns
midas_filtered <- midas_taxonomy %>% 
  select(kingdom, phylum, class, order, family, genus, species, species_id)

# Join in the taxonomy columns into metaG
e0026_metaG_full <- metaGdata %>% 
  left_join(midas_filtered, by="species_id")

# Save full dataset into a text file
write.table(e0026_metaG_full,
            file = paste0(OUTDIR,"e0026_metaG_full.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)



