# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-getFullyLostFams/100825-getFullyLostFams.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-testNullModelLostFamilies/out/"

COMPONENTS <- c("0", "1", "2", "3")
PALETTE.COMPONENTS <- c(  "cornflowerblue","coral2", "darkorange","chartreuse3")
names(PALETTE.COMPONENTS) <- COMPONENTS


curr_replicate <- 1


# --------------------------- Try out one community (XBB + XBA) ----------------------------
# 
# # Get day 36 mix otus
# XBAXHB_day36_mixture <- actual_colonizers_results %>% 
#   filter(day == "036", subject1 == "XBA", subject2 == "XHB") %>% 
#   select(-colonized_day29, colonized_day36, colonized_day64)

XBB_donor_day36 <- single_donor_ASVs %>% 
  filter(biosample1 == "XBB-029")

# Get day 29 colonizers
XBAXHB_day29_colonizers <- actual_colonizers_results %>% 
  filter(day == "029", subject2 == "XHB", actual_colonizer == 1)

# Get day 36 colonizers
XBAXHB_day36_colonizers <- actual_colonizers_results %>% 
  filter(day == "036", subject2 == "XHB", actual_colonizer == 1)

# Get day 36 recipient otus
XBA_recipient <- recipient_ASVs %>% 
  filter(biosample1 == "XBA-036")

# Get lost otus
XBA_lost_strains <- fully_lost_families %>% 
  filter(subject1 == "XBA")

# Add a component classification for each otu
XBAXHB_day36_mixture_components <- XBAXHB_day36_mixture %>% 
  mutate(component = case_when(
    OTU %in% XBA_recipient$OTU ~ "1", # OTUs present in post-abx recipient
    OTU %in% XBAXHB_day29_colonizers$OTU ~ "2", # OTUs present as colonizers in day 29 AND day 36 mixtures
    Family %in% XBA_lost_strains$Family ~ "3", # OTUs lost post-abx in the recipient
    TRUE ~ "0" # "other"
  ))


