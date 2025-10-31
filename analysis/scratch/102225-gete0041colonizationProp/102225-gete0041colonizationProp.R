# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102625-gete0041LostTaxa/102625gete0041LostTaxa.R")


get_colonization <- function(mix_id, recipient_id, donor_id) {
  mix_otus <- e0041.A %>% 
    filter(mixture == mix_id, replicate == curr_replicate, recipient != "blank", donor != "blank")
  
  recipient_otus <- e0041_control_recipients %>% 
    filter(recipient ==  recipient_id) # replicate specified in loading script. Currently replicate == 1
  
  donor_otus <- e0041_control_donors %>%
    filter(donor == donor_id) # replicate specified in loading script. Currently replicate == 1
  
  # potential colonizers are those present in the donor but not in the recipient
  potential_colonizers <- donor_otus$OTU[!(donor_otus$OTU %in% recipient_otus$OTU)]
  
  donor_potential_otus <- donor_otus %>% 
    filter(OTU %in% potential_colonizers)
  
  # colonizers are the subset of potentials that did make it into the mix
  mixture_colonization <- mix_otus %>% 
    mutate(colonizer = ifelse(OTU %in% potential_colonizers, 1, 0))
  
  # return both the dataframe and potential colonizer list
  return(list(
    mixture_colonization = mixture_colonization,
    donor_potential_otus = donor_potential_otus
  ))
}

# mix_id <- "pre-abx+XBB-029"
# recipient_id <- "pre-abx"
# donor_id <- "XBB-029"

# run on all mixtures: 

mixture_colonization_full <- tibble()
potential_colonizers_full <- tibble()
mixture_colonization_proportion_full <- tibble()

# Loop over mixtures
for (mix_id in mixture_ids) {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[2]
  recipient_id <- ids[1]
  
  # run colonization function w/ potential colonizer dataframe too
  result <- get_colonization(mix_id, recipient_id, donor_id) # get colonization for a specific mixture
  
  # extract each of the dataframes
  mixture_colonization <- result$mixture_colonization # get colonizers
  potential_colonizers <- result$donor_potential_otus # get potential colonizers
  
  # count the number of potential and actual colonizers
  n_potential_colonizers <- nrow(potential_colonizers)
  n_colonizers <- mixture_colonization %>%
    filter(colonizer == 1) %>% 
    nrow()
  
  # get values for a single mixture
  mixture_row <- tibble(
    mixture = mix_id,
    recipient = recipient_id,
    donor = donor_id,
    n_potential_colonizers = n_potential_colonizers,
    n_colonizers = n_colonizers,
    colonization_proportion = n_colonizers / n_potential_colonizers
  )
  
  # bind all the mixture rows together
  mixture_colonization_full <- bind_rows(mixture_colonization_full, mixture_colonization)
  potential_colonizers_full <- bind_rows(potential_colonizers_full, potential_colonizers)
  mixture_colonization_proportion_full <- bind_rows(mixture_colonization_proportion_full, mixture_row)
}

# ---------------------------------------- ADD DIFFERENTIAL COLONIZER LABEL

# We need to group OTUs by their mixture pair, independent of day
diff_colonizers <- mixture_colonization_full %>%
  group_by(OTU) %>% 
  summarize(
    colonized_pre_abx = as.integer(any(recipient == "pre-abx" & colonizer == 1)),
    colonized_post_abx_v1 = as.integer(any(recipient == "post-abx-V1" & colonizer == 1)),
    colonized_post_abx_v2 = as.integer(any(recipient == "post-abx-V2" & colonizer == 1)),
    .groups = "drop"
  ) %>%
  mutate(diff_colonizer_v1 = as.integer((colonized_post_abx_v1) & !colonized_pre_abx), diff_colonizer_v2 = as.integer((colonized_post_abx_v2) & !colonized_pre_abx))


mixture_colonization_full <- mixture_colonization_full %>%
  left_join(diff_colonizers,
            by = c("OTU"))

# --------------------------------------- ADD LOST FAMILY LABEL
# mixtures
mixture_colonization_full <- mixture_colonization_full %>% 
  mutate(lost_V1 = ifelse(Family %in% post_abx_lost_families_v1, 1, 0), lost_V2 = ifelse(Family %in% post_abx_lost_families_v2, 1, 0))

# recipients
e0041_control_recipients <- e0041_control_recipients %>% 
  mutate(lost_V1 = ifelse(Family %in% post_abx_lost_families_v1, 1, 0), lost_V2 = ifelse(Family %in% post_abx_lost_families_v2, 1, 0))


# --------------------------------------- FIND WHO "OTHER" COLONIZERS ARE, MATCH AT THE ORDER LEVEL?

other_colonizers_v1 <- mixture_colonization_full %>% 
  filter(recipient == "post-abx-V1", diff_colonizer_v1 ==1, lost_V1 ==0) %>%
  group_by(mixture, Family) %>% 
  summarize(family_counts = n())

lost_colonizers_v1 <- unique(mixture_colonization_full %>% 
  filter(recipient == "post-abx-V1", diff_colonizer_v1 ==1, lost_V1 ==1))


# --------------------------------------- WHO DOESN'T COLONIZE?

mixture_colonization_pre_abx <- mixture_colonization_full %>% 
  filter(recipient == "pre-abx")
mixture_colonization_post_abx_V1 <- mixture_colonization_full %>% 
  filter(recipient == "post-abx-V1")
mixture_colonization_post_abx_V2 <- mixture_colonization_full %>% 
  filter(recipient == "post-abx-V2")

e0041_control_donors <- e0041_control_donors %>% 
  mutate(colonized_pre_abx = ifelse(OTU %in% mixture_colonization_pre_abx$OTU, 1, 0), 
         colonized_post_abx_V1 = ifelse(OTU %in% mixture_colonization_post_abx_V1$OTU, 1, 0),
         colonized_post_abx_V2 = ifelse(OTU %in% mixture_colonization_post_abx_V2$OTU, 1, 0))

# Investigate colonizers
colonizers_V1 <- e0041_control_donors %>% 
  filter(colonized_post_abx_V1 == 1) %>% 
  group_by(donor, Family) %>% 
  summarize(family_counts = n())

colonizers_V2 <- e0041_control_donors %>% 
  filter(colonized_post_abx_V2 == 1) %>% 
  group_by(donor, Family) %>% 
  summarize(family_counts = n())

# Investigate interesting features about taxa that don't colonize
e0041_control_recipients_post_abx_v1 <- e0041_control_recipients %>% 
  filter(recipient == "post-abx-V1")

not_colonizers_V1 <- e0041_control_donors %>% 
  filter(colonized_post_abx_V1 == 0) %>% 
  mutate(shared_recipient = ifelse(OTU %in% e0041_control_recipients_post_abx_v1$OTU, 1, 0)) %>% 
  filter(shared_recipient == 0) %>% 
  group_by(donor, Family) %>% 
  summarize(family_counts = n())

e0041_control_recipients_pre_abx <- e0041_control_recipients %>% 
  filter(recipient == "pre-abx")

not_colonizers_V1 <- e0041_control_donors %>% 
  filter(colonized_pre_abx == 0) %>% 
  mutate(shared_recipient = ifelse(Family %in% e0041_control_recipients_pre_abx$Family, 1, 0)) %>% 
  filter(shared_recipient == 0) %>% 
  group_by(donor, Family) %>% 
  summarize(family_counts = n())






