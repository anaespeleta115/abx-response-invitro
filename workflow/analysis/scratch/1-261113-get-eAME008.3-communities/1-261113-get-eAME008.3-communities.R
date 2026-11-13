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


# Get paths to all CSV files in the 4YRV6C folder
files <- list.files(
  path = "~/Documents/GitHub/abx-response-invitro/data/4YRV6C_results",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Read and combine into one dataframe
bottomUp_data <- files |>
  set_names(basename(files)) |>
  map_dfr(read_csv, .id = "sample") |> 
  mutate(community = str_sub(sample, -22, -20))

# Adjust formatting, filter out replicate 2 communities
bottomUp_data <- bottomUp_data |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family)) |> 
  filter(!community %in% c("B01", "B02", "B03"))



