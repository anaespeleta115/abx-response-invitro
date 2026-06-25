# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
library(foreach)


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/out/"


# ------------------------ COMPUTE COLONIZATION EFFICACY ---------------------------

get_colonization <- function(mix_id, donor_id, recipient_id, replicate) {
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == curr_replicate, relAbundance > limit_of_detection)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == curr_replicate, relAbundance > limit_of_detection)
  
  donor_asvs_subset <- single_donor_ASVs %>%  
    filter(biosample == donor_id, relAbundance > limit_of_detection)
  
  # potential colonizers are those present in the donor but not in the recipient
  potential_colonizers <- donor_asvs_subset$OTU[!(donor_asvs_subset$OTU %in% recipient_asvs_subset$OTU)]
  
  donor_asvs_potential <- donor_asvs_subset %>% 
    filter(OTU %in% potential_colonizers)
  
  # colonizers are the subset of potentials that actually did make it into the mix
  mix_asvs_colonization <- mix_asvs_subset %>% 
    mutate(colonizer = ifelse(OTU %in% potential_colonizers, 1, 0))
  
  # return both the dataframe and potential colonizer list
  return(list(
    colonization_df = mix_asvs_colonization,
    donor_asvs_potential = donor_asvs_potential
  ))
}



# --------------- LOOP OVER ALL MIXTURES AND MAKE FINAL COLONIZATION DATAFRAMES ----------------

colonization_results <- tibble()
potential_colonizer_results <- tibble()
colonization_prop_results <- tibble()


# Loop over mixtures
for (mix_id in mixture_ids) {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[1]
  recipient_id <- ids[2]
  
  result <- get_colonization(mix_id, donor_id, recipient_id, 1) # get colonization for a specific mixture
  
  colonizers <- result$colonization_df # get colonizers
  potential_colonizers <- result$donor_asvs_potential # get potential colonizers
  
  n_potential_colonizers <- nrow(potential_colonizers)
  n_colonizers <- colonizers %>%
    filter(colonizer == 1) %>% 
    nrow()
  
  # bind all the mixture rows together
  colonization_results <- bind_rows(colonization_results, colonizers)
  potential_colonizer_results <- bind_rows(potential_colonizer_results, potential_colonizers)
  
  prop_row <- tibble(
    mixture = mix_id,
    biosample2 = donor_id,
    biosample1 = recipient_id,
    n_potential_colonizers = n_potential_colonizers,
    n_colonizers = n_colonizers,
    prop_colonizers = n_colonizers / n_potential_colonizers
  )
  
  # bind all the mixture rows together
  colonization_prop_results <- bind_rows(colonization_prop_results, prop_row)
  
}

# add additional columns for the colonization proportion dataset
colonization_prop_results <- colonization_prop_results %>%
  mutate(
    recipient = str_sub(biosample1, 1, -5),
    donor = str_sub(biosample2, 1, -5),
    day = str_sub(biosample1, -3)
  )


# ---------------- GET DIFFERENTIAL COLONIZERS ----------------------

colonization_results <- colonization_results %>% 
  mutate(mixture_pair = str_replace_all(mixture, "-\\d{3}", ""))

# Group OTUs by their mixture pair, independent of their day
diff_colonizers <- colonization_results %>%
  group_by(mixture_pair, OTU) %>%
  summarize(
    colonized_day29 = as.integer(any(day == "029" & colonizer == 1)),
    colonized_day36 = as.integer(any(day == "036" & colonizer == 1)),
    colonized_day64 = as.integer(any(day == "064" & colonizer == 1)),
    .groups = "drop"
  ) %>%
  mutate(diff_colonizer_36 = as.integer((colonized_day36) & !colonized_day29), diff_colonizer_64 = as.integer((colonized_day64) & !colonized_day29 & !colonized_day36))  # Should we also consider differential colonizers that don't invade day 36 but do by day 64


colonization_results <- colonization_results %>%
  left_join(diff_colonizers,
            by = c("mixture_pair", "OTU"))


# Write out base dataset with differential colonizer flags into a text file

write_delim(colonization_results, paste0(OUTDIR, "differential_colonization_AME004.txt"))
