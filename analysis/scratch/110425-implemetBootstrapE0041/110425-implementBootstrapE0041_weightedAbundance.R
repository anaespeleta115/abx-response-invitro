# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-gete0041colonization/102125-gete0041colonization.R")


# Set seed to keep the same randomization
set.seed(123)

# ------------------------------------- Get post-abx V1 mixtures -------------------------------------------
mixture_colonization_post_abx1 <- mixture_colonization_full %>%
  filter(recipient == "pre-abx" | recipient == "post-abx-V1") 


# ------------------------------------- Get post-abx V2 mixtures -------------------------------------------
mixture_colonization_post_abx2 <- mixture_colonization_full %>%
  filter(recipient == "pre-abx" | recipient == "post-abx-V2") 




# ------------------------------------- Define a function to get potential colonizers

get_potential_colonizers <- function(mix, donor_id, recipient_id) {
  
  mix_subset <- mixture_colonization_full %>% 
    filter(mixture == mix)
  
  recipient_subset <- e0041_control_recipients %>% 
    filter(recipient == recipient_id)
  
  donor_subset <- e0041_control_donors %>%  
    filter(donor == donor_id)
  
  # potential colonizers are those present in the donor but not in the recipient
  potential_colonizers <- donor_subset$OTU[!(donor_subset$OTU %in% recipient_subset$OTU)]
  
  donor_potential <- donor_subset %>% 
    filter(OTU %in% potential_colonizers)
  
  # return potential colonizer list
  return(donor_potential)
}


# --------------------------- Bootstrap setup --------------------------------

# Initialize results list
bootstrap_results_list <- list()

mixture_ids <- mixture_colonization_post_abx1 %>%
  filter(donor != "super", recipient == "post-abx-V1") %>% 
  distinct(mixture) %>%
  pull(mixture)

# Loop over each mixture
for (mix in mixture_ids) {
  
  ids <- unlist(strsplit(mix, "\\+"))
  donor_id <- ids[2]
  recipient_id <- ids[1]
  
  # Get all colonizers for that mixture
  sample_colonizers <- mixture_colonization_post_abx1 %>%
    filter(mixture == mix, recipient == "post-abx-V1", colonized_post_abx_v1 == 1) %>%
    distinct(OTU, relAbundance, Family)
  
  # Get all differential colonizers for that mixture
  sample_diff_colonizers <- mixture_colonization_post_abx1 %>%
    filter(donor == donor_id, diff_colonizer_v1 == 1) %>%
    distinct(OTU, relAbundance, Family)
  
  # Get lost taxa between pre-abx and post-abx-V1
  recipient_lost <- e0041_control_recipients %>%
    filter(recipient == "pre-abx", lost_V1 == 1) %>% 
    distinct(OTU, relAbundance, Family)
  
  # Get all potential colonizers for that mixture
  sample_potential_colonizers <- get_potential_colonizers(mix, donor_id, recipient_id)
  
  # --------------------------- Run Bootstrap --------------------------------
  n_trials <- 1000
  n_sample <- nrow(sample_diff_colonizers) # I think we want the pool of colonizers to be the size of differential colonizers so that both
  # the observed and bootstrap pools are the same size.
  
  bootstrap_results <- map_dfr(1:n_trials, function(trial) {
    boot_sample <- sample_colonizers %>%
      slice_sample(n = n_sample, weight_by = relAbundance, replace = TRUE) # take a slice out of the post-abx-v1 colonizers and see if they
  # are enriched for lost taxa
    
    otu_shared <- length(intersect(boot_sample$OTU, recipient_lost$OTU))
    fam_shared_weighted <- sum(boot_sample$relAbundance[boot_sample$Family %in% recipient_lost$Family])
    
    otu_shared_ids <- intersect(boot_sample$OTU, recipient_lost$OTU)
    # fam_shared_ids <- intersect(boot_sample$Family, recipient_lost$Family)
    
    # Replace empty vectors with "none"
    otu_shared_ids <- if (length(otu_shared_ids) == 0) "none" else otu_shared_ids
    # fam_shared_ids <- if (length(fam_shared_ids) == 0) "none" else fam_shared_ids
    
    tibble(trial = trial, 
           shared_otus = otu_shared, 
           shared_families = fam_shared_weighted, 
           otu_ids = list(otu_shared_ids))
  })
  
  # Observed values: raw count of OTUs/Families shared between real diff colonizers and real lost taxa
  observed_otus <- length(intersect(sample_diff_colonizers$OTU, recipient_lost$OTU))
  observed_families <- sum(sample_diff_colonizers$relAbundance[sample_diff_colonizers$Family %in% recipient_lost$Family])
  
  observed_otu_ids <- intersect(sample_diff_colonizers$OTU, recipient_lost$OTU)
  # observed_families_ids <- intersect(sample_diff_colonizers$Family, recipient_lost$Family)
  
  # Replace empty vectors with "none"
  observed_otu_ids <- if (length(observed_otu_ids) == 0) "none" else observed_otu_ids
  # observed_families_ids <- if (length(observed_families_ids) == 0) "none" else observed_families_ids
  
  # Add observed values as a final row
  bootstrap_results <- bootstrap_results %>%
    mutate(trial = as.character(trial)) %>%
    bind_rows(tibble(
      trial = "Observed",
      shared_otus = observed_otus,
      shared_families = observed_families,
      otu_ids = list(observed_otu_ids)
      # fam_ids = list(observed_families_ids)
    ))
  
  # Store for this mixture
  bootstrap_results_list[[mix]] <- bootstrap_results
}

combined_bootstrap_results <- bind_rows(bootstrap_results_list, .id = "mixture")

# Extract observed rows per mixture
observed_rows <- combined_bootstrap_results %>%
  filter(trial == "Observed") %>%
  select(mixture, observed_otus = shared_otus, observed_families = shared_families,
         observed_otu_ids = otu_ids)

# Join observed values to bootstrap results
enrichment_pvals <- combined_bootstrap_results %>%
  filter(trial != "Observed") %>%
  left_join(observed_rows, by = "mixture") %>%
  group_by(mixture) %>%
  summarise(
    p_value_fam = mean(shared_families >= observed_families),
    p_value_otu = mean(shared_otus >= observed_otus)
  )

