source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")


# Compute species richness

donor_richness <- single_donor_ASVs %>% 
  filter(relAbundance > limit_of_detection) %>% 
  group_by(replicate, passage, biosample1) %>% 
  summarize(richness = n())