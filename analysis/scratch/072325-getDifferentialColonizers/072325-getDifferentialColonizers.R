# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072225-getColonizationProportion/072225-getColonizationProportion.R")


actual_colonizers_results <- actual_colonizers_results %>% 
  mutate(day = str_sub(biosample1, -3), 
         mixture_pair = str_replace_all(mixture, "-\\d{3}", ""))

# We need to group OTUs by their mixture pair, independent of day
diff_colonizers <- actual_colonizers_results %>%
  group_by(mixture_pair, OTU) %>%
  summarize(
    colonized_day29 = as.integer(any(day == "029" & actual_colonizer == 1)),
    colonized_day36 = as.integer(any(day == "036" & actual_colonizer == 1)),
    colonized_day64 = as.integer(any(day == "064" & actual_colonizer == 1)),
    .groups = "drop"
  ) %>%
  mutate(diff_colonizer_36 = as.integer((colonized_day36) & !colonized_day29), diff_colonizer_64 = as.integer((colonized_day64) & !colonized_day29))  # Should we also consider differential colonizers that don't invade day 36 but do by day 64


actual_colonizers_results <- actual_colonizers_results %>%
  left_join(diff_colonizers %>% select(mixture_pair, OTU, diff_colonizer_36, diff_colonizer_64),
            by = c("mixture_pair", "OTU"))