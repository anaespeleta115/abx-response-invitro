# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")

curr_replicate <- 1

# Compute colonization as an efficacy metric and return a colonization dataframe as well as the potential colonizer list

get_colonization <- function(mix_id, donor_id, recipient_id, replicate) {
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
colonization_prop_results <- tibble()


# Loop over mixtures
for (mix_id in mixture_ids) {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[1]
  recipient_id <- ids[2]
  
  result <- get_colonization(mix_id, donor_id, recipient_id, 1) # get colonization for a specific mixture
  
  actual_colonizers <- result$colonization_df # get colonizers
  potential_colonizers <- result$donor_asvs_potential # get potential colonizers
  
  n_potential_colonizers <- nrow(potential_colonizers)
  n_actual_colonizers <- actual_colonizers %>%
    filter(actual_colonizer == 1) %>% 
    nrow()
  
  # bind all the mixture rows together
  actual_colonizers_results <- bind_rows(actual_colonizers_results, actual_colonizers)
  potential_colonizers_results <- bind_rows(potential_colonizers_results, potential_colonizers)
  
  prop_row <- tibble(
    mixture = mix_id,
    biosample2 = donor_id,
    biosample1 = recipient_id,
    n_potential_colonizers = n_potential_colonizers,
    n_actual_colonizers = n_actual_colonizers,
    prop_colonizers = n_actual_colonizers / n_potential_colonizers
  )
  
  # bind all the mixture rows together
  colonization_prop_results <- bind_rows(colonization_prop_results, prop_row)
  
}


colonization_prop_results <- colonization_prop_results %>%
  mutate(
    recipient = str_sub(biosample1, 1, -5),
    donor = str_sub(biosample2, 1, -5),
    day = str_sub(biosample1, -3),
    household = str_sub(biosample1, 1, -6)
  )


