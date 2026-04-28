source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")



# Get recipient communities
recipient_ASVs <- eAME004_data %>% 
  filter(relAbundance > limit_of_detection) %>%    
  filter(biosample2 == "blank") %>% 
  filter(biosample1 != "blank" )


# Extract donor community IDs from e0029
donor_communities <- eAME004_data %>%
  filter(biosample2 != "blank") %>% 
  pull(biosample2) %>%
  unique()

# Extract single donor OTUs from e0026
# ADD BACK FILTER: relAbundance > limit_of_detection
single_donor_ASVs <- eACE010_data %>%
  filter(relAbundance > limit_of_detection,
         biosample %in% donor_communities)


# Extract subject, day, household information and filter out extra donor tests
eAME004_base <- eAME004_data %>%
  filter(
    biosample1 != "blank",
    biosample2 != "blank",
    biosample2 != "B-mix",
    !str_detect(biosample2, "\\+"))


# Get mixture communities
mixture_ASVs <- eAME004_base %>% 
  filter(relAbundance > limit_of_detection)

# Get list of mixture IDs and  out b-mix and two-mix donors
mixture_ids <- unique(mixture_ASVs %>%
                        pull(mixture))

