library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

# set global variables
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-261110-get-eAME008-compositions/out/"

# get data from Plasmidsaurus output (taxonomy table) for community in well A01
A01_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_1_AME0108_A01/3W4XCX_1_AME0108_A01_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "1")


# plot composition of well A01, plate eAME0108, community 1, replicate 1.
p_A01_community_composition <- A01_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 1, Rep 1", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_A01_community_composition

savePNGPDF(paste0(OUTDIR, "p_A01_community_composition"), p_A01_community_composition, 3, 3)



# get data from Plasmidsaurus output (taxonomy table) for community in well B01
B01_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_4_AME0108_B01/3W4XCX_4_AME0108_B01_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "1")


# plot composition of well B01, plate eAME0108, community 1, replicate 2.
p_B01_community_composition <- B01_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 1, Rep 2", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_B01_community_composition

savePNGPDF(paste0(OUTDIR, "p_B01_community_composition"), p_B01_community_composition, 3, 3)



# get data from Plasmidsaurus output (taxonomy table) for community in well A02
A02_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_2_AME0108_A02/3W4XCX_2_AME0108_A02_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "2")


# plot composition of well A02, plate eAME0108, community 2, replicate 1.
p_A02_community_composition <- A02_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 2, Rep 1", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_A02_community_composition

savePNGPDF(paste0(OUTDIR, "p_A02_community_composition"), p_A02_community_composition, 3, 3)




# get data from Plasmidsaurus output (taxonomy table) for community in well B02
B02_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_5_AME0108_B02/3W4XCX_5_AME0108_B02_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "2")


# plot composition of well B02, plate eAME0108, community 2, replicate 2.
p_B02_community_composition <- B02_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 2, Rep 2", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_B02_community_composition

savePNGPDF(paste0(OUTDIR, "p_B02_community_composition"), p_B02_community_composition, 3, 3)



# get data from Plasmidsaurus output (taxonomy table) for community in well A03
A03_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_3_AME0108_A03/3W4XCX_3_AME0108_A03_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "3")


# plot composition of well A03, plate eAME0108, community 3, replicate 1.
p_A03_community_composition <- A03_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 3, Rep 1", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_A03_community_composition

savePNGPDF(paste0(OUTDIR, "p_A03_community_composition"), p_A03_community_composition, 3, 3)



# get data from Plasmidsaurus output (taxonomy table) for community in well B03
B03_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_6_AME0108_B03/3W4XCX_6_AME0108_B03_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "3")


# plot composition of well B03, plate eAME0108, community 3, replicate 2.
p_B03_community_composition <- B03_community |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Comm 3, Rep 2", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_B03_community_composition

savePNGPDF(paste0(OUTDIR, "p_B03_community_composition"), p_B03_community_composition, 3, 3)