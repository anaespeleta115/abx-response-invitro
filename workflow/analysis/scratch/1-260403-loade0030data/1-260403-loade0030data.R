library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library("phyloseq")
library(PNWColors)
library(paletteer)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260403-loade0030data/out/"

household_data <- read.table("~/Documents/GitHub/abx-response-invitro/data/e0026-e0029-e0030.txt", header = TRUE)

# Separate household data into each experiment
e0030 <- household_data %>% 
  filter(experiment == "e0030")

e0029 <- household_data %>% 
  filter(experiment == "e0029")

# Define global variables
curr_replicate <- 1

limit_of_detection <-  1e-3

my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 
