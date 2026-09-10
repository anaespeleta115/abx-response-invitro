library(tidyverse)
library(dplyr)
library(ggplot2)
library("phyloseq")
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

# Define global variables

limit_of_detection <-  1e-3

my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-2608012-loade0018data/out/"

e0018_data <- fread("~/Documents/GitHub/abx-response-invitro/data/e0018_ps_all.txt.gz", header = TRUE)
e0018_ps <- readRDS("~/Documents/GitHub/abx-response-invitro/data/e0018_ps_all.rds")

# COMMENT/UNCOMMENT here to switch sample sheets
# e0018_sample_sheet <- read_csv("~/Documents/GitHub/abx-response-invitro/data/e0018_ss_updated_20220223.csv") %>%
#   mutate(household = ifelse(mixture == "C to A-pre" | mixture == "C to A-post", str_sub(biosample2, 2, 2), household))
e0018_sample_sheet <- read_csv("~/Documents/GitHub/abx-response-invitro/data/e0018_ss_updated_20210901.csv")

# Define global variables

limit_of_detection <-  1e-3

my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


# join in metadata from sample sheet
e0018_data_full <- e0018_data %>% 
  left_join(e0018_sample_sheet %>% select(sample, mixture, household, biosample_full), by = "sample")


# Get only autologous mixtures

e0018_autologous_mixes <- e0018_data_full %>% 
  filter(mixture == "A-pre" | mixture == "A-post" | mixture == "A-pre to A-post") %>% 
  mutate(mix_type = paste(mixture, inoculationReplicate, sep = "-")) %>% 
  filter(!(mixture == "A-pre to A-post" & ratio == "1:1")) %>% 
  mutate(mix_type = fct_relevel(mix_type, "A-pre-1", "A-post-1", "A-pre to A-post-1", "A-pre to A-post-2", "A-pre to A-post-3"))


# Get pre-abx recipient OTUs
preabx_recipient_OTUs <- e0018_autologous_mixes %>% 
  filter(relAbundance > limit_of_detection) %>%    
  filter(mixture == "A-pre")

# Get post-abx recipient OTUs
postabx_recipient_OTUs <- e0018_autologous_mixes %>% 
  filter(relAbundance > limit_of_detection) %>%    
  filter(mixture == "A-post")

# Get mixture OTUs
mixture_OTUs <- e0018_autologous_mixes %>% 
  filter(relAbundance > limit_of_detection) %>%    
  filter(mixture == "A-pre to A-post" & mix_type == "A-pre to A-post-1") # only use one replicate for now

# Get non-redundant list of mixture IDs
mixture_ids <- unique(mixture_OTUs %>%
                        pull(biosample_full))


# Get only cohabiting mixtures + recipients and donors

e0018_cohabiting_mixes <- e0018_data_full %>% 
  filter(mixture == "A-pre" | mixture == "A-post" | mixture == "B" | mixture == "B to A-pre" | mixture == "B to A-post") %>% 
  mutate(mix_type = paste(mixture, inoculationReplicate, sep = "-")) %>% 
  mutate(mix_type = fct_relevel(mix_type, "A-pre-1", "A-post-1", "B to A-pre-1", "B to A-pre-2", "B to A-pre-3", 
                                "B to A-post-1", "B to A-post-2", "B to A-post-3"))

# Get only unrelated mixtures + recipients and donors

e0018_unrelated_mixes <- e0018_data_full %>% 
  filter(mixture == "A-pre" | mixture == "A-post" | mixture == "C to A-pre" | mixture == "C to A-post") %>% 
  mutate(mix_type = paste(mixture, inoculationReplicate, sep = "-")) %>% 
  mutate(mix_type = fct_relevel(mix_type, "A-pre-1", "A-post-1", "C to A-pre-1", "C to A-pre-2", "C to A-pre-3", 
                                "C to A-post-1", "C to A-post-2", "C to A-post-3"))

