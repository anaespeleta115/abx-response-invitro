# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-plotSpeciesRichnessMG/031726-plotSpeciesRichnessMG.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-plotSpeciesRichnessHousehold/031726-plotSpeciesRichnessHousehold.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/workflow/analysis/plotDefaults.R")
# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-getMGSpeciesChange/out"







# Change day variables to possible column names
MG_p8_species_richness <- MG_p8_species_richness %>% 
  mutate(day = case_when(
    str_detect(day, "001") ~ "day_1",
    str_detect(day, "029") ~ "day_29",
    str_detect(day, "036") ~ "day_36",
    str_detect(day, "064") ~ "day_64",
    TRUE ~ "0"
  ))

# Find the change in species richness across all subjects, filter for antibiotic-taking subjects
MG_delta_diversity <- MG_p8_species_richness %>% 
  pivot_wider(names_from = day, values_from = species_richness) %>% 
  group_by(subject) %>% 
  mutate(delta_1_64 = day_64 - day_1, delta_29_36 = day_36 - day_29) %>% 
  mutate(loss_1_64 = ifelse(delta_1_64 < 0, 1, 0), loss_29_36 = ifelse(delta_29_36 < 0, 1, 0)) %>% 
  filter(subject %in% subjectsAbx)


write.csv(MG_delta_diversity, "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-getMGSpeciesChange/out/MG_delta_diversity", row.names = FALSE)







# Change day variables to possible column names
MG_p0_species_richness <- MG_p0_species_richness %>% 
  mutate(day = case_when(
    str_detect(day, "001") ~ "day_1",
    str_detect(day, "029") ~ "day_29",
    str_detect(day, "036") ~ "day_36",
    str_detect(day, "064") ~ "day_64",
    TRUE ~ "0"
  ))

# Find the change in species richness across all subjects, filter for antibiotic-taking subjects
household_delta_diversity <- MG_p0_species_richness %>% 
  pivot_wider(names_from = day, values_from = species_richness) %>% 
  group_by(subject) %>% 
  mutate(delta_1_64 = day_64 - day_1, delta_29_36 = day_36 - day_29) %>% 
  mutate(loss_1_64 = ifelse(delta_1_64 < 0, 1, 0), loss_29_36 = ifelse(delta_29_36 < 0, 1, 0)) %>% 
  filter(subject %in% subjectsAbx)


write.csv(MG_delta_diversity, "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-getMGSpeciesChange/out/MG_delta_diversity.csv", row.names = FALSE)


write.csv(household_delta_diversity, "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-getMGSpeciesChange/out/household_data_diverstiy.csv", row.names = FALSE)

