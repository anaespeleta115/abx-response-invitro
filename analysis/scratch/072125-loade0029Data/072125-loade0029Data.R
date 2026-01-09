# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library(purrr)
library(patchwork)
library(foreach)
library(forcats)
library(cowplot)
library(stringr)
library(scico)
library(paletteer)


# Set global variables
my_colors <- readRDS("/Users/aespelet/Documents/Github/abx-response-invitro/data/familyColorPalette.rds") 

# path reference: /Users/aespelet/Documents/Github/

# Set limit of detection
limit_of_detection <- 1e-3

# Read in household dataset
household_data <- read.table("/Users/aespelet/Documents/Github/abx-response-invitro/data/e0026-e0029-e0030.txt", header = TRUE)

# Divide dataset into separate tables by experiment
e0029 <- filter(household_data, str_detect(household_data$sample, "e0029"))
e0026 <- filter(household_data, str_detect(household_data$sample, "e0026"))


# Filter out duplicated controls
e0029 <- e0029 %>%
  filter(
    !(biosample2 == "blank") |
      (biosample2 == "blank" & !well %in% c("G11", "H11"))) %>%
  filter(
    biosample1 != "blank" |
      (plate == "e0029.B.5" & well %in% c("G12", "H12"))
  )

# Get recipient, donor, and mixture communities into separate datasets

# Adjust replicate!!!
replicate <- 1

# Get recipient communities
recipient_ASVs <- e0029 %>% 
  filter(relAbundance > limit_of_detection, replicate == replicate) %>%    
  filter(biosample2 == "blank") %>% 
  filter(
    !(biosample2 == "blank") |
      (biosample2 == "blank" & !well %in% c("G11", "H11"))) %>%
  filter(biosample1 != "blank" )


# Extract donor community IDs from e0029
donor_communities <- e0029 %>%
  pull(biosample2) %>%
  unique()

# Extract single donor OTUs from e0026
# ADD BACK FILTER: relAbundance > limit_of_detection
single_donor_ASVs <- e0026 %>%
  filter(biosample1 %in% donor_communities, passage == "8")


# Extract subject, day, household information and filter out extra donor tests
e0029_base <- e0029 %>%
  filter(
    biosample1 != "blank",
    biosample2 != "blank",
    biosample2 != "B-mix",
    !str_detect(biosample2, "\\+")
  ) %>%
  mutate(
    subject1 = str_sub(biosample1, 1, -5),
    subject2 = str_sub(biosample2, 1, -5),
    day = str_sub(biosample1, -3),
    household = str_sub(biosample1, 1, -6)
  )


# Get mixture communities
mixture_ASVs <- e0029_base %>% 
  filter(relAbundance > limit_of_detection, replicate == replicate)


# Get list of mixture IDs and  out b-mix and two-mix donors
mixture_ids <- unique(mixture_ASVs %>%
                        pull(mixture))


mixture_ids_36 <- unique(mixture_ASVs %>% filter(day == "036") %>% 
                        pull(mixture)) 

# Set palettes

PASSAGE <- c(0, 8)
PALETTE.PASSAGE <- c(  "gray80", paletteer_d("nationalparkcolors::Saguaro"))
names(PALETTE.PASSAGE) <- PASSAGE

SUBJECT <- c("XBA", "XDA", "XEA", "XKA")
PALETTE.SUBJECT <- c(paletteer_d("nationalparkcolors::Badlands"))
names(PALETTE.SUBJECT) <- SUBJECT

DAY <- c("0", "7", "35")
PALETTE.DAY <- c(paletteer_d("nationalparkcolors::Voyageurs"))
names(PALETTE.DAY) <- DAY




