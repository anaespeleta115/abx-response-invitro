# load scripts 
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102025-loade0040data/102025-loade0040data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-gete0040-e0041diversity/out/"


e0040_diversity <- e0040 %>% 
  filter(replicate == 1, relAbundance > limit_of_detection) %>% 
  group_by(community) %>% 
  summarize(species_richness = n())


e0041_diversity <- e0041.A %>% 
  filter(replicate == 1, relAbundance > limit_of_detection) %>% 
  group_by(mixture) %>% 
  summarize(species_richness = n())

