library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

# set global variables
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 

limit_of_detection <- 1e-3

# get data from Plasmidsaurus output (taxonomy table) for community in well A01
A01_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_1_AME0108_A01/3W4XCX_1_AME0108_A01_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "1") 

# get data from Plasmidsaurus output (taxonomy table) for community in well B01
B01_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_4_AME0108_B01/3W4XCX_4_AME0108_B01_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "1")

# ---------------------------------


# get data from Plasmidsaurus output (taxonomy table) for community in well A02
A02_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_2_AME0108_A02/3W4XCX_2_AME0108_A02_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "2")


# get data from Plasmidsaurus output (taxonomy table) for community in well B02
B02_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_5_AME0108_B02/3W4XCX_5_AME0108_B02_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "2")

# ---------------------------------

# get data from Plasmidsaurus output (taxonomy table) for community in well A03
A03_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_3_AME0108_A03/3W4XCX_3_AME0108_A03_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "3")

# get data from Plasmidsaurus output (taxonomy table) for community in well B03
B03_community <- read_csv("~/Documents/GitHub/abx-response-invitro/data/3W4XCX_results/3W4XCX_6_AME0108_B03/3W4XCX_6_AME0108_B03_taxonomy-table.csv") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  mutate(community = "3")





