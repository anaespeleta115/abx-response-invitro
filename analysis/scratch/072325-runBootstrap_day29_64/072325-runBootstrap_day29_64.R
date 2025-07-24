# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")


set.seed(123)

# Exclude day 64 data
actual_colonizers_results_filtered <- actual_colonizers_results %>%
  mutate(subject = str_sub(biosample1, 1, -5)) %>% 
  filter(day %in% c("029", "064")) %>%
  mutate(day = str_sub(biosample1, -3)) 

# Initialize results list
bootstrap_results_list <- list()

mixture_ids <- actual_colonizers_results_filtered %>%
  distinct(mixture) %>%
  pull(mixture)

# Loop over each mixture
for (mix in mixture_ids) {
  
  ids <- unlist(strsplit(mix, "\\+"))
  donor_id <- ids[1]
  donor_id <- str_sub(donor_id, 1, -5)
  recipient_id <- ids[2] 
  recipient_id <- str_sub(recipient_id, 1, -5)
  mix_pair <- paste(donor_id, recipient_id, sep = "+")
  
  
  # Get all colonizers for that mixture
  single_sample_actual_colonizers <- actual_colonizers_results_filtered %>%
    filter(mixture == mix, day %in% c("064"), actual_colonizer == 1) %>%
    distinct(OTU, Family)
  
  # Get all differential colonizers for that mixture
  single_sample_diff_colonizers <- actual_colonizers_results_filtered %>%
    filter(mixture_pair == mix_pair, diff_colonizer_64 == 1) %>%
    distinct(OTU, Family)
  
  # Get lost taxa between day 29 and day 64
  single_recipient_lost <- recipient_lost_29_64 %>%
    mutate(subject = str_sub(biosample1, 1, -5)) %>% 
    filter(subject == recipient_id) %>%
    distinct(OTU, Family)
  
  # If not enough data, skip
  if (nrow(single_sample_actual_colonizers) == 0 || nrow(single_recipient_lost) == 0) {
    next
  }
  
  # Bootstrap setup
  n_trials <- 1000
  n_sample <- nrow(single_sample_actual_colonizers)
  
  bootstrap_results <- map_dfr(1:n_trials, function(trial) {
    boot_sample <- single_sample_actual_colonizers %>%
      slice_sample(n = n_sample, replace = TRUE)
    
    otu_shared <- length(intersect(boot_sample$OTU, single_recipient_lost$OTU))
    fam_shared <- length(intersect(boot_sample$Family, single_recipient_lost$Family))
    
    otu_shared_ids <- intersect(boot_sample$OTU, single_recipient_lost$OTU)
    fam_shared_ids <- intersect(boot_sample$Family, single_recipient_lost$Family)
    
    # Replace empty vectors with "none"
    otu_shared_ids <- if (length(otu_shared_ids) == 0) "none" else otu_shared_ids
    fam_shared_ids <- if (length(fam_shared_ids) == 0) "none" else fam_shared_ids
    
    tibble(trial = trial, shared_otus = otu_shared, shared_families = fam_shared, otu_ids = list(otu_shared_ids) , fam_ids = list(fam_shared_ids))
  })
  
  # Observed values: raw count of OTUs/Families shared between real diff colonizers and real lost taxa
  observed_otus <- length(intersect(single_sample_diff_colonizers$OTU, single_recipient_lost$OTU))
  observed_families <- length(intersect(single_sample_diff_colonizers$Family, single_recipient_lost$Family))
  
  observed_otu_ids <- intersect(single_sample_diff_colonizers$OTU, single_recipient_lost$OTU)
  observed_families_ids <- intersect(single_sample_diff_colonizers$Family, single_recipient_lost$Family)
  
  # Replace empty vectors with "none"
  observed_otu_ids <- if (length(observed_otu_ids) == 0) "none" else observed_otu_ids
  observed_families_ids <- if (length(observed_families_ids) == 0) "none" else observed_families_ids
  
  # Add observed values as a final row
  bootstrap_results <- bootstrap_results %>%
    mutate(trial = as.character(trial)) %>%
    bind_rows(tibble(
      trial = "Observed",
      shared_otus = observed_otus,
      shared_families = observed_families,
      otu_ids = list(observed_otu_ids),
      fam_ids = list(observed_families_ids)
    ))
  
  # Store for this mixture
  bootstrap_results_list[[mix]] <- bootstrap_results
}

combined_bootstrap_results <- bind_rows(bootstrap_results_list, .id = "mixture")


# Extract observed rows per mixture
observed_rows <- combined_bootstrap_results %>%
  filter(trial == "Observed") %>%
  select(mixture, observed_otus = shared_otus, observed_families = shared_families,
         observed_otu_ids = otu_ids, observed_family_ids = fam_ids)

# Calculate enrichment p-values from bootstrap rows
enrichment_pvals <- combined_bootstrap_results %>%
  filter(trial != "Observed") %>%
  mutate(shared_families = as.numeric(shared_families)) %>%
  group_by(mixture) %>%
  summarise(
    p_value_fam = mean(shared_families >= observed_rows$observed_families[observed_rows$mixture == unique(mixture)]),
    p_value_otu = mean(shared_otus >= observed_rows$observed_otus[observed_rows$mixture == unique(mixture)]),
    .groups = "drop"
  )

# Merge observed info with p-values
enrichment_summary_64 <- left_join(observed_rows, enrichment_pvals, by = "mixture") %>%
  mutate(
    enriched_fam = p_value_fam < 0.05,
    enriched_otu = p_value_otu < 0.05
  )



ids_64 <- enrichment_summary_64 %>%
  mutate(
    observed_otu_ids = sapply(observed_otu_ids, paste, collapse = ", "),
    observed_family_ids = sapply(observed_family_ids, paste, collapse = ", ")
  ) %>%
  select(mixture, observed_otu_ids, observed_family_ids)



