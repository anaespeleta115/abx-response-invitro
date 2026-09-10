source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260817-calculate-e0018-quality/out/"


# -----------
## This script calculates for each autologous mixture in e0018 the percentage of OTUs that are present in either the pre- or post-
## antibiotic recipient
# -----------


# ----------- here I calculate the percentage for one household ---------------

A_pre_abx <- preabx_recipient_OTUs %>% 
  filter(household == "A")

A_post_abx <- postabx_recipient_OTUs %>% 
  filter(household == "A")

A_mixtures <- mixture_OTUs %>% 
  filter(household == "A")

A_e0018_purity <- A_mixtures %>% 
  mutate(impure = ifelse(!(OTU %in% A_pre_abx$OTU |
                           OTU %in% A_post_abx$OTU), 
                         1, 0)) %>%
  mutate(impure_otus = sum(impure), total_otus = n()) %>% 
  group_by(household, total_otus, impure_otus) %>% 
  summarize(impurity_percentage = sum(impure)/n()*100)
  


# -------- here I calculate this percentage across all households ------------------------


# get a vector of all unique households present in mixture_OTUs
households <- unique(mixture_OTUs$household)

# initialize an empty data frame to store results
all_households_purity <- data.frame()

# loop through each household
for (hh in households) {
  
  # filter datasets for the current household iteration
  curr_pre_abx  <- preabx_recipient_OTUs  %>% filter(household == hh)
  curr_post_abx <- postabx_recipient_OTUs %>% filter(household == hh)
  curr_mixtures <- mixture_OTUs           %>% filter(household == hh)
  
  # calculate purity metrics for the current household
  curr_summary <- curr_mixtures %>%
    mutate(
      impure = if_else(
        !(OTU %in% curr_pre_abx$OTU | OTU %in% curr_post_abx$OTU),
        1, 0
      )
    ) %>%
    summarize(
      household           = hh,
      total_otus          = n(),
      impure_otus         = sum(impure),
      pure_otus           = total_otus - impure_otus,
      impurity_percentage = (impure_otus / total_otus) * 100,
      purity_percentage   = 100 - impurity_percentage
    )
  
  # append current household result to the master dataset
  all_households_purity <- bind_rows(all_households_purity, curr_summary)
}

# mutate a column to specify source samplesheet
all_households_purity <- all_households_purity %>% 
  mutate(samplesheet = "20210901")

# export dataset as a csv file
write.csv(all_households_purity, "household_purity_20210901.csv", row.names = FALSE)

