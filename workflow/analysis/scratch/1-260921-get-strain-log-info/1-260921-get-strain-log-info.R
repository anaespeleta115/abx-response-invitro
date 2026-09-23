library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(data.table)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

# set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260921-get-strain-log-info/out/"

# set global variables
limit_of_detection <-  1e-3
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 


# read in strain log
strain_log <- read_csv("~/Documents/GitHub/abx-response-invitro/data/Xue_lab_strainLog_v2 - strainLog.csv") |> 
  filter(!is.na(TAX_totalReads), !is.na(TAX_species...4), !TAX_family...3 == "Pyronemataceae", BIO_subject == "XBA" |BIO_subject == "XBA") |> # filter out unnecessary info
  rename(TAX_purity = TAX_purity...5, TAX_family = TAX_family...3, 
         TAX_species = TAX_species...4) |> # rename columns
  select(-TAX_family...25, -TAX_species...27, -TAX_purity...29, -nickname) |> # remove repeated columns
  mutate(TAX_family = ifelse(TAX_species == "Phocaeicola vulgatus", "Bacteroidaceae", TAX_family),
         TAX_family = ifelse(TAX_species == "Flintibacter sp. KGMB00164", "Oscillospiraceae", TAX_family),
         TAX_family = ifelse(TAX_species == "Flavonifractor plautii", "Oscillospiraceae", TAX_family))  # assign correct family to unassigned species


# list out column names
colnames(strain_log)


# try plotting number of isolates per family
family_distribution <- strain_log |> 
  group_by(TAX_family) |> 
  summarize(unique_species_per_fam = n_distinct(TAX_species), total_species_per_fam = n())


p_family_distribution <- family_distribution |> 
  ggplot(aes(x = fct_reorder(TAX_family, -unique_species_per_fam), y = unique_species_per_fam, fill = TAX_family)) +
  geom_col(color = "black") +
    labs(
      x = "Family",
      y = "Num. total isolates",
      fill = "Family"
    ) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 10)) +
    DEFAULTS.THEME_PRINT +
    theme(
      axis.text.x = element_text(
        hjust = 1, vjust = 0.5, size = 6, angle = 90
      ),
      legend.position = "none"
    )

p_family_distribution

savePNGPDF(paste0(OUTDIR, "p_family_distribution_unique"), p_family_distribution, 3, 3)


# count number of species in total
num_species <- strain_log |> 
  summarize(total_species = n_distinct(TAX_species), total_families = n_distinct(TAX_family))
  

# subset log to keep only relevant info
unique_strains <- strain_log |> 
  select(strainID, TAX_family, TAX_genus, TAX_species, TAX_purity, ISO_plate, LOC_culture, LOC_glycerol) |> 
  group_by(TAX_species) |>
  slice_max(order_by = TAX_purity, n = 1, with_ties = FALSE) # with ties = FALSE ensures that we only keep 1 even if they have the same value.


# # Save dataset into a text file
# write.table(unique_strains ,
#             file = paste0(OUTDIR,"unique_strains.txt"),
#             sep = "\t",
#             row.names = FALSE,
#             quote = FALSE)



# plot the purity distribution
purity_density <- unique_strains |> 
  ggplot(aes(x = TAX_purity)) +
  geom_density() +
  xlim(0, NA)+
  labs(
    title = "",
    x = "Purity",
    y = "Density"
  ) +
  DEFAULTS.THEME_PRINT

purity_density

savePNGPDF(paste0(OUTDIR, "purity_density_unique"), purity_density , 3, 3)


p_num_media_types <- num_media_types |> 
  ggplot(aes(x = fct_reorder(ISO_plate, -num_species), y = num_species, fill = ISO_plate)) +
  geom_col(color = "black") +
  labs(
    x = "Media type",
    y = "Num. total isolates"
  ) +
  scale_y_continuous(limits = c(0, 30)) +
  DEFAULTS.THEME_PRINT +
  theme(
    axis.text.x = element_text(
      hjust = 1, vjust = 0.5, size = 6, angle = 90
    ),
    legend.position = "none"
  )

p_num_media_types

savePNGPDF(paste0(OUTDIR, "p_num_media_types_unique"), p_num_media_types , 3, 3)




num_plates <- unique_strains |> 
  mutate(LOC_glycerol_plateID = str_sub(LOC_glycerol, 1, 7)) |> 
  group_by(LOC_glycerol_plateID) |> 
  summarize(num_species_per_plate = n())

p_num_plates <- num_plates |> 
  ggplot(aes(x = fct_reorder(LOC_glycerol_plateID, -num_species_per_plate), y = num_species_per_plate, fill = LOC_glycerol_plateID)) +
  geom_col(color = "black") +
  labs(
    x = "Glycerol Plate ID",
    y = "Num. total isolates"
  ) +
  scale_y_continuous(limits = c(0, 15)) +
  DEFAULTS.THEME_PRINT +
  theme(
    axis.text.x = element_text(
      hjust = 1, vjust = 0.5, size = 6, angle = 90
    ),
    legend.position = "none"
  )

p_num_plates

savePNGPDF(paste0(OUTDIR, "p_num_plates_unique"), p_num_plates , 3, 3)




  