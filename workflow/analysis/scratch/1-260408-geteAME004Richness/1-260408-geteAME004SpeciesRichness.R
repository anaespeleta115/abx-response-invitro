source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-geteAME004Richness/out/"


# Leave only mixtures with "blank", recipients by themselves.
eAME004_data_recipients <- eAME004_data %>% 
  filter(biosample2 == "blank")
  
  
# Calculate the diversity in each well.
eAME004_recipient_diversity <- eAME004_data_recipients %>%
  group_by(subject, day) %>%
  filter(relAbundance>0.001) %>%
  summarize(nASVs=n())


# Calculate change in species richness

# Change day variables to possible column names
eAME004_recipient_diversity <- eAME004_recipient_diversity %>% 
  mutate(day = case_when(
    str_detect(day, "029") ~ "day_29",
    str_detect(day, "036") ~ "day_36",
    str_detect(day, "064") ~ "day_64",
    TRUE ~ "0"
  ))

# Find the change in species richness across all subjects, filter for antibiotic-taking subjects
eAME004_delta_diversity <- eAME004_recipient_diversity %>% 
  pivot_wider(names_from = day, values_from = nASVs) %>% 
  group_by(subject) %>% 
  mutate(delta_29_36 = day_36 - day_29) %>% 
  mutate(loss_29_36 = ifelse(delta_29_36 < 0, 1, 0))



# eACE010 pooled communities
write.csv(eAME004_delta_diversity,
            file = paste0(OUTDIR,"eACE010_full_pooled.csv"),
            row.names = FALSE,
            quote = FALSE)



# 
# 
# # Plot a stacked bar plot of composition in a select set of communities.
# eAME004_recipient_plot <- eAME004_data_recipients %>%
#   filter(biosample1 %in% c("XDA-029", "XDA-036", "XDA-064"),
#          relAbundance>0.001) %>%
#   ggplot() +
#   geom_bar(aes(x=mixture, y=relAbundance, fill=Family), stat="identity", color="black")
# 
# eAME004_recipient_plot