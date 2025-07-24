# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")


get_colonization <- function(mix_id, donor_id, recipient_id, replicate) {
  
  # Subset mixture and recipient ASVs for the given replicate
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == replicate)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == replicate)
  
  donor_asvs_subset <- single_donor_ASVs %>%  # change this when we have donor_asvs with mixed donors
    filter(biosample1 == donor_id)
  
  # Get shared OTUs between donor and mixture.  Add colonization column: 1 if OTU is shared with donor, 0 otherwise
  mix_asvs_colonization <- mix_asvs_subset %>% 
    mutate(colonization = ifelse(OTU %in% donor_asvs_subset$OTU & !(OTU %in% recipient_asvs_subset$OTU), 1, 0))
  
  return(mix_asvs_colonization)
}


# Apply get_colonization() function across all mixture IDs
colonization_results <- foreach(mix_id = mixture_ids, .combine = bind_rows) %do% {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[1]
  recipient_id <- ids[2]
  
  result <- get_colonization(
    mix_id = mix_id,
    donor_id = donor_id,
    recipient_id = recipient_id,
    replicate = 1
  )
  
  result %>% mutate(mixture = mix_id)
}

