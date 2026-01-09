source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102025-loade0040data/102025-loade0040data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")



# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102625-gete0041LostTaxa/out/"


e0040_pre_abx <- unique(e0040 %>% 
  filter(replicate == 1, community == "XEA-pre-abx") %>% 
  pull(OTU))

e0040_postabx_v1 <- unique(e0040 %>% 
  filter(replicate == 1, community == "XEA-post-abx-V1") %>% 
  pull(OTU))

e0040_postabx_v2 <- unique(e0040 %>% 
  filter(replicate == 1, community == "XEA-post-abx-V2") %>% 
  pull(OTU))

# post_abx_lost_taxa_v1 <- setdiff(e0040_pre_abx, e0040_postabx_v1)
# post_abx_lost_taxa_v2 <- setdiff(e0040_pre_abx, e0040_postabx_v2)


e0041_pre_abx <- unique(e0041_control_recipients %>% 
                          filter(recipient == "pre-abx") %>% 
                          pull(OTU))

e0041_postabx_v1 <- unique(e0041_control_recipients %>% 
                             filter(recipient == "post-abx-V1") %>% 
                             pull(OTU))

e0041_postabx_v2 <- unique(e0041_control_recipients %>% 
                             filter(recipient == "post-abx-V2") %>% 
                             pull(OTU))

post_abx_lost_taxa_v1 <- setdiff(e0040_pre_abx, e0040_postabx_v1)
post_abx_lost_taxa_v2 <- setdiff(e0040_pre_abx, e0040_postabx_v2)

e0041_recipients_lost_V1 <- e0041_control_recipients %>% filter(OTU %in% post_abx_lost_taxa_v1) %>%
  select(OTU, Family) %>%
  distinct()

e0041_recipients_lost_V2 <- e0041_control_recipients %>% filter(OTU %in% post_abx_lost_taxa_v2) %>%
  select(OTU, Family) %>%
  distinct()