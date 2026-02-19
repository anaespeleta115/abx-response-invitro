# Load data

source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/110725-getE0026SpeciesChange/out/"


# Set palettes
pal <- pnw_palette("Sunset2", 2, type = "discrete")

ABX <- c("No", "Yes")
PALETTE.ABX <- c("gray80","#88CCEE")
names(PALETTE.ABX) <- ABX


PASSAGE <- c(0, 8)
PALETTE.PASSAGE <- c(  "gray80", "#B084CC")
names(PALETTE.PASSAGE) <- PASSAGE


# Get species richness for all samples in passage 0
e0026_diversity_p0 <- e0026 %>% 
  filter(relAbundance > limit_of_detection, passage == 0) %>% 
  group_by(passage, subject, day) %>% 
  summarize(species_richness = n())


# Get species richness for all samples in passage 8
e0026_diversity_p8 <- e0026 %>% 
  filter(relAbundance > limit_of_detection, passage == 8) %>% 
  group_by(passage, subject, day) %>% 
  summarize(species_richness = n())

# Get species richness for all samples
e0026_diversity_all <- e0026 %>% 
  filter(relAbundance > limit_of_detection) %>% 
  group_by(passage, subject, day) %>% 
  summarize(species_richness = n())


# Change day variables to possible column names
e0026_diversity_p8 <- e0026_diversity_p8 %>% 
  mutate(day = case_when(
    str_detect(day, "001") ~ "day_1",
    str_detect(day, "029") ~ "day_29",
    str_detect(day, "036") ~ "day_36",
    str_detect(day, "064") ~ "day_64",
    TRUE ~ "0"
  ))

# Find the change in species richness across all subject, filter for antibiotic-taking subjects
e0026_delta_diversity <- e0026_diversity_p8 %>% 
  pivot_wider(names_from = day, values_from = species_richness) %>% 
  group_by(subject) %>% 
  mutate(delta_1_64 = day_64 - day_1, delta_29_36 = day_36 - day_29) %>% 
  mutate(loss_1_64 = ifelse(delta_1_64 < 0, 1, 0), loss_29_36 = ifelse(delta_29_36 < 0, 1, 0)) %>% 
  filter(subject %in% subjectsAbx)


write.csv(e0026_delta_diversity, "~/Documents/GitHub/abx-response-invitro/analysis/scratch/110725-getE0026SpeciesChange/out/e0026_delta_diversity.csv", row.names = FALSE)
