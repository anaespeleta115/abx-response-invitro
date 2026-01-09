source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/111425-loade0058data/111425-loade0058data.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/111623-verifyCompositions/out/"


# Figure out what was actually dropped between two communities
full_community <- e0058 %>% 
  filter(str_detect(community, "full_community"), version == "all_strains") %>% 
  select(OTU, Family)

no_bacter_community <- e0058 %>% 
  filter(community == "bacteroides_dropped_2") %>%
  select(OTU, Family)

dropped <- setdiff(full_community, no_bacter_community) # get those strains in full community, but not
# in the bacteroides_dropped community


# dropped Bifido
bifido_dropped_communities <- e0058 %>% 
  filter(str_detect(community, "bifido")) %>% 
  group_by(community) %>% 
  filter(Family == "Bifidobacteriaceae") %>% 
  summarize(num_bifido_strains = n_distinct(OTU))

