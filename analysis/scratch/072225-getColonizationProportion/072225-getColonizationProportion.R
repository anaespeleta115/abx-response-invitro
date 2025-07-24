# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")



# Compute colonization as an efficacy metric and return a colonization dataframe as well as the potential colonizer list

get_colonization <- function(mix_id, donor_id, recipient_id, replicate) {
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == replicate)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == replicate)
  
  donor_asvs_subset <- single_donor_ASVs %>%  
    filter(biosample1 == donor_id)
  
  # potential colonizers are those present in the donor but not in the recipient
  potential_colonizers <- donor_asvs_subset$OTU[!(donor_asvs_subset$OTU %in% recipient_asvs_subset$OTU)]
  
  donor_asvs_potential <- donor_asvs_subset %>% 
    filter(OTU %in% potential_colonizers)
  
  # actual colonizers are the subset of potentials that did make it into the mix
  mix_asvs_colonization <- mix_asvs_subset %>% 
    mutate(actual_colonizer = ifelse(OTU %in% potential_colonizers, 1, 0))
  
  # return both the dataframe and potential colonizer list
  return(list(
    colonization_df = mix_asvs_colonization,
    donor_asvs_potential = donor_asvs_potential
  ))
}



# Flag potential colonizers that do colonize as "actual colonizers"

actual_colonizers_results <- tibble()
potential_colonizers_results <- tibble()

# Loop over mixtures
for (mix_id in mixture_ids) {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[1]
  recipient_id <- ids[2]
  
  result <- get_colonization(mix_id, donor_id, recipient_id, 1)
  
  
  actual_colonizers <- result$colonization_df 
  potential_colonizers <- result$donor_asvs_potential
  
  
  # bind all the mixture rows together
  actual_colonizers_results <- bind_rows(actual_colonizers_results, actual_colonizers)
  potential_colonizers_results <- bind_rows(potential_colonizers_results, potential_colonizers)
}

