library(data.table)
library(tidyverse)
library(ggplot2)
library("phyloseq")
library(cowplot)
library(patchwork)

source("~/Documents/GitHub/abx-response-invitro/analysis/plotDefaults.R")

e0058 <- fread("~/Documents/GitHub/abx-response-invitro/data/e0058ps_all.txt.gz")

e0058_obj <- readRDS("~/Documents/GitHub/abx-response-invitro/data/ps_all_e0058.rds")


# Define color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 

# Define replicate used for colonization function
curr_replicate <- 1

# set limit of detection
limit_of_detection <- 0.001

# get a recipient column
e0058 <- e0058 %>% 
  mutate(
    community = case_when(
      str_detect(well, "A1$|E1$") ~ "full_community",
      str_detect(well, "B1$|F1$") ~ "full_community_ss",
      str_detect(well, "A2$|E2$") ~ "bacteroides_dropped",
      str_detect(well, "B2$|F2$") ~ "bacteroides_dropped_ss",
      str_detect(well, "A3|E3") ~ "bifido_dropped",
      str_detect(well, "B3|F3") ~ "bifido_dropped_ss",
      str_detect(well, "A4|E4") ~ "hydro_dropped",
      str_detect(well, "B4|F4") ~ "hydro_dropped_ss",
      str_detect(well, "A5|E5") ~ "lachno_dropped",
      str_detect(well, "B5|F5") ~ "lachno_dropped_ss",
      str_detect(well, "A6|E6") ~ "prevo_dropped",
      str_detect(well, "B6|F6") ~ "prevo_dropped_ss",
      str_detect(well, "A7|E7") ~ "rumino_dropped",
      str_detect(well, "B7|F7") ~ "rumino_dropped_ss",
      str_detect(well, "A8|E8") ~ "strep_dropped",
      str_detect(well, "B8|F8") ~ "strep_dropped_ss",
      str_detect(well, "A9|A10$|E9|E10$") ~ "blank",
      str_detect(well, "B9|B10$|F9|F10$") ~ "blank",
      str_detect(well, "A11$|E11$") ~ "strep_dropped_half_fam_1",
      str_detect(well, "B11$|F11$") ~ "strep_dropped_half_fam_2",
      str_detect(well, "A12$|E12$") ~ "bacteroides_dropped_half_fam_1",
      str_detect(well, "B12$|F12$") ~ "bacteroides_dropped_half_fam_2",
      
      str_detect(well, "C1$|G1$") ~ "bacteroides_dropped_1",
      str_detect(well, "D1$|H1$") ~ "lachno_dropped_2",
      str_detect(well, "C2$|G2$") ~ "bacteroides_dropped_2",
      str_detect(well, "D2$|H2$") ~ "lachno_dropped_3",
      str_detect(well, "C3|G3") ~ "bacteroides_dropped_3",
      str_detect(well, "D3|H3") ~ "peptS_dropped_1",
      str_detect(well, "C4|G4") ~ "bacteroides_dropped_4",
      str_detect(well, "D4|H4") ~ "rumino_dropped_1",
      str_detect(well, "C5|G5") ~ "bifido_dropped_1",
      str_detect(well, "D5|H5") ~ "rumino_dropped_2",
      str_detect(well, "C6|G6") ~ "bifido_dropped_2",
      str_detect(well, "D6|H6") ~ "strep_dropped_1",
      str_detect(well, "C7|G7") ~ "bifido_dropped_3",
      str_detect(well, "D7|H7") ~ "strep_dropped_2",
      str_detect(well, "C8|G8") ~ "bifido_dropped_4",
      str_detect(well, "D8|H8") ~ "strep_dropped_3",
      str_detect(well, "C9|G9") ~ "erysipelotrichaceae_dropped_1",
      str_detect(well, "D9|H9") ~ "strep_dropped_4",
      str_detect(well, "C10$|G10$") ~ "erysipelotrichaceae_dropped_2",
      str_detect(well, "D10$|H10$") ~ "strep_dropped_5",
      str_detect(well, "C11$|G11$") ~ "familyXI_dropped_1",
      str_detect(well, "C12$|G12$") ~ "lachno_dropped_1",
      str_detect(well, "D11$|H11$|D12$|H12$") ~ "blank",
      TRUE ~ NA_character_
    ))

# add more specific version column

e0058 <- e0058 %>% 
  mutate(
    replicate = case_when(
      str_detect(well, "A|B") ~ "1",
      str_detect(well, "E|F") ~ "2",
      str_detect(well, "C|D") ~ "1",
      str_detect(well, "D|H") ~ "2",
      
      # str_detect(well, "A11$|E11$") ~ "half_fam_1",
      # str_detect(well, "B11$|F11$") ~ "half_fam_2",
      # str_detect(well, "A12$|E12$") ~ "half_fam_1",
      # str_detect(well, "B12$|F12$") ~ "half_fam_2",
      TRUE ~ NA_character_
    ))


# get ordered wells
e0058 <- e0058 %>%
  mutate(well = sprintf("%s%02d",
                        substr(well, 1, 1),
                        as.integer(substring(well, 2))))

rows <- LETTERS[1:8]
cols <- sprintf("%02d", 1:12)
wells <- as.vector(sapply(rows, function(r) paste0(r, cols)))

# filter for the current replicate
e0058 <- e0058 %>% 
  filter(replicate == 1)
