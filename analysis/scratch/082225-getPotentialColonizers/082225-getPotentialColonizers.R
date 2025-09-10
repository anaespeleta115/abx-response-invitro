# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")



curr_replicate <- 1

# mix_id <- "XBB-029+XDA-029"
# 
# recipient_id <- "XDA-029"
# 
# donor_id <- "XBB-029"


get_potential_colonizers <- function(mix_id, donor_id, recipient_id, replicate) {
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == curr_replicate)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == curr_replicate)
  
  donor_asvs_subset <- single_donor_ASVs %>%  
    filter(biosample1 == donor_id)
  
  # potential colonizers are those present in the donor but not in the recipient
  potential_colonizers <- donor_asvs_subset$OTU[!(donor_asvs_subset$OTU %in% recipient_asvs_subset$OTU)]
  
  donor_asvs_potential <- donor_asvs_subset %>% 
    filter(OTU %in% potential_colonizers)
  
  # return potential colonizer list
  return(donor_asvs_potential)
}

# poten <- get_potential_colonizers(mix_id = mix_id, recipient_id = recipient_id, donor_id = donor_id, replicate = curr_replicate)

