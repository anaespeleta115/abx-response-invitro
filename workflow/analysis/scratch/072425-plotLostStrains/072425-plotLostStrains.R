# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")




# Get lost strains between days 29 and 36

num__fam_lost_36 <- recipient_lost_29_36  %>% 
  distinct(Family, biosample1) %>%  # remove duplicate Family–biosample1 combos
  count(Family, name = "n_samples")


# Get lost strains between days 36 and 64

num_fam_lost_64 <- recipient_lost_29_64  %>% 
  distinct(Family, biosample1) %>%  # remove duplicate Family–biosample1 combos
  count(Family, name = "n_samples")


# Number of strains lost per subject day 36

num_lost_strains_29_36 <- recipient_lost_29_36 %>% 
  group_by(biosample1) %>% 
  summarize(n_OTU = n_distinct(OTU))


  
# Number of strains lost per subject day 64

num_lost_strains_29_64 <- recipient_lost_29_64 %>% 
  group_by(biosample1) %>% 
  summarize(n_OTU = n_distinct(OTU))


