source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260903-percent-recovered-e0018/out/"


## This script calculates the percent of OTUs recovered in the autologous mixture out of all lost ones in the A-post recipient

### First, we have to calculate all lost OTUs in the recipient


# # extract day 29 taxa not in day 36 taxa
# lost_OTUs <- preabx_recipient_OTUs |>  
#   group_by(household) |> 
#   filter(!(OTU %in% postabx_recipient_OTUs$OTU)) |>   # get all OTUs not present in the post-abx community
#   mutate(lost_OTU = 1)




# get a vector of all unique households present in mixture_OTUs
households <- unique(mixture_OTUs$household)

# initialize an empty data frame to store results
all_households_recovered_OTUs <- data.frame()

# loop through each household
for (hh in households) {
  
  # filter datasets for the current household iteration
  curr_pre_abx  <- preabx_recipient_OTUs  %>% filter(household == hh)
  curr_post_abx <- postabx_recipient_OTUs %>% filter(household == hh)
  curr_mix <- mixture_OTUs %>% filter(household == hh)
  
  # add lost flag to pre-abx recipient OTUs
  curr_pre_abx_lost <- curr_pre_abx %>%
    mutate(
      lost = if_else(
        !(OTU %in% curr_post_abx$OTU),
        1, 0
      )
    )
  
  # add gained flag to pre-abx recipient OTUs 
  curr_pre_abx_recovered <- curr_pre_abx_lost |> 
    mutate(
      recovered = if_else(
        (lost == 1 & OTU %in% curr_mix$OTU),
        1, 0
      )
    ) %>%
    summarize(
      household           = hh,
      lost_otus          = sum(lost),
      recovered_otus           = sum(recovered),
      recovery_percentage = (recovered_otus / lost_otus) * 100
    )
  
  
  
  
  # append current household result to the master dataset
  all_households_recovered_OTUs <- bind_rows(all_households_recovered_OTUs, curr_pre_abx_recovered)
}

# export dataset as a csv file
write.csv(all_households_recovered_OTUs, "all_households_recovered_OTUs.csv", row.names = FALSE)



