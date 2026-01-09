source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/111425-loade0058data/111425-loade0058data.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/111623-verifyCompositions/out/"


# full communities species richness
full_community_richness <- e0058 %>% 
  filter(str_detect(community, "full_community"), relAbundance > limit_of_detection) %>% 
  group_by(community) %>% 
  summarize(species_richness = n())

full_community_richness_by_fam <- e0058 %>% 
  filter(str_detect(community, "full_community"), relAbundance > limit_of_detection) %>% 
  group_by(community, Family) %>% 
  summarize(species_richness_fam = n())

# dropped Bacteroides
bacter_dropped_communities <- e0058 %>% 
  filter(str_detect(community, "bacter")) %>% 
  group_by(community) %>% 
  filter(Family == "Bacteroidaceae") %>% 
  summarize(num_bacter_strains = n_distinct(OTU))

# Figure out what was actually dropped between two communities -------------------------

full_community <- e0058 %>% 
  filter(str_detect(community, "full_community")) %>% 
  select(OTU, Family)

# with 
no_bacter_community <- e0058 %>% 
  filter(community == "bacteroides_dropped") %>%
  select(OTU, Family)

dropped_all <- setdiff(full_community, no_bacter_community) # get those strains in full community but not
# in the bacteroides_dropped community


no_bacter_community <- e0058 %>% 
  filter(community == "bacteroides_dropped_ss") %>%
  select(OTU, Family)

dropped_all_ss

