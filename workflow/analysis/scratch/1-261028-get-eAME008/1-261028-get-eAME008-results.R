library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-261028-get-eAME008-results/out/"

# set global variables
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


eAME008_results <- read.table("~/Documents/GitHub/abx-response-invitro/data/eAME008_top-hits.txt", header = TRUE)

taxonomy_table <- read.table("~/Documents/GitHub/abx-response-invitro/data/eAME008_taxonomy.txt", header = TRUE)
