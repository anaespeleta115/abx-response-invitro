# load libraries

library(data.table)
library(tidyverse)
library(ggplot2)
library("phyloseq")
library(cowplot)
library(patchwork)
library(paletteer)

# scale_colour_paletteer_d("lisa::BridgetRiley")
# scale_color_paletteer_d("lisa::BridgetRiley")
# scale_fill_paletteer_d("lisa::BridgetRiley")
# paletteer_d("lisa::BridgetRiley")

# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/plotDefaults.R")

e0040 <- fread("~/Documents/GitHub/bottom-up-community-mixtures/data/ps_all.txt.gz")

e0040_obj <- readRDS("~/Documents/GitHub/bottom-up-community-mixtures/data/ps_all.rds")

# Define color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 

# set limit of detection
limit_of_detection <- 1e-3

e0040_RCM <- e0040 %>% 
  filter(well %in% c("E2","F2","G2","H2"))

# assign community label and replicate based on wells
e0040 <- e0040 %>% 
  mutate(replicate = case_when(
    well %in% c("A1","A2","A3") ~ "1", 
    well %in% c("B1","B2","B3") ~ "2",
    well %in% c("C1","C2","C3") ~ "3", 
    well %in% c("D1","D2","D3") ~ "4"
  )) %>% 
  mutate(community = case_when(
    well %in% c("A1","B1","C1", "D1") ~ "XEA-pre-abx", 
    well %in% c("A2","B2","C2", "D2") ~ "XEA-post-abx-V1",
    well %in% c("A3","B3","C3", "D3") ~ "XEA-post-abx-V2",
  )) %>% 
  filter(!is.na(community))


e0040_control_recipients_richness <- e0040 %>% group_by(replicate, community, well) %>% summarize(richness = sum(relAbundance>0.001))

