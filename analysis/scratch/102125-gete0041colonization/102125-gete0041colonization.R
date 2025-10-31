# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
library(foreach)

get_colonization <- function(mixture_id, recipient_id, donor_id) {
  
  # Subset mixture and recipient ASVs for the given replicate
  mix_otus <- e0041.A %>% 
    filter(mixture == mixture_id, replicate == curr_replicate, recipient != "blank", donor != "blank")
  
  recipient_otus <- e0041_control_recipients %>% 
    filter(recipient ==  recipient_id) # replicate specified in loading script. Currently replicate == 1

  donor_otus <- e0041_control_donors %>%
    filter(donor == donor_id) # replicate specified in loading script. Currently replicate == 2
  
  # Get shared OTUs between donor and mixture.  Add colonization column: 1 if OTU is shared with donor, 0 otherwise
  mixture_colonization <- mix_otus %>% 
    mutate(
      colonizer = ifelse(OTU %in% donor_otus$OTU & !(OTU %in% recipient_otus$OTU), 1, 0))

  return(mixture_colonization)
}

pre_abx_colonization <- get_colonization("pre-abx+XBB-029", "pre-abx", "XBB-029")


# Apply get_colonization() function across all mixture IDs
XEA_colonization <- foreach(mix_id = mixture_ids, .combine = bind_rows) %do% {
  ids <- unlist(strsplit(mix_id, "\\+"))
  donor_id <- ids[1]
  recipient_id <- ids[2]
  
  result <- get_colonization(
    mix_id,
    recipient_id,
    donor_id
  )
}