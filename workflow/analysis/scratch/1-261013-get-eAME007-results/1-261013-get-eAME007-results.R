library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-261013-get-eAME007-results/out/"

# set global variables
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


eAME007_results <- read_csv("~/Documents/GitHub/abx-response-invitro/data/eAME007-annotated-results.csv") |> 
  mutate(well_padded = str_sub(`plate-well`, 9, 11))

taxonomy_table <- read.table("~/Documents/GitHub/abx-response-invitro/data/eAME007-2_taxonomy.txt", header = TRUE)

# join in taxonomy information and eAME007 results
eAME007_full_results <- eAME007_results |> 
  left_join(taxonomy_table, by = "well_padded") |> 
  mutate(Family = ifelse(Family == "Unknown", "Oscillospiraceae", Family))


# try plotting number of isolates per family
family_distribution <- eAME007_full_results |> 
  filter(purity == 1, !is.na(Species)) |> # filter for strains with purity of 1
  group_by(Family) |> 
  summarize(unique_species_per_fam = n_distinct(Species), total_species_per_fam = n())

# plot families by the number of species
p_family_distribution <- family_distribution |> 
  ggplot(aes(x = fct_reorder(Family, -unique_species_per_fam), y = unique_species_per_fam, fill = Family)) +
  geom_col(color = "black") +
  labs(
    x = "Family",
    y = "Num. total isolates",
    fill = "Family"
  ) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 5)) +
  DEFAULTS.THEME_PRINT +
  theme(
    axis.text.x = element_text(
      hjust = 1, vjust = 0.5, size = 6, angle = 90
    ),
    legend.position = "none"
  )

p_family_distribution

savePNGPDF(paste0(OUTDIR, "p_family_distribution_unique"), p_family_distribution, 3, 3)


unique_isolates <- family_distribution |> 
  summarize(total_isolates = sum(unique_species_per_fam))


