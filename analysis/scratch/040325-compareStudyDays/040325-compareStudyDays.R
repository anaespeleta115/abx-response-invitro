# Plot stacked bar plots with a facet for each subject, showing compositional changes over time across timepoints, at passage 8.

# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/040325-compareStudyDays/out/"


### Compositional changes pre and post-abx

P8_compositions <- combined_day_data %>% 
  filter(passage == 8)

# Summarize relative abundance per biosample1 by Family
composition_pre_postAbx <- P8_compositions %>%
  group_by(biosample1, Family) %>%
  summarise(relAbundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    P8_compositions %>% distinct(biosample1, day, subject, household, antibiotic),
    by = "biosample1"
  ) 
# 
# %>%
#   mutate(
#     day = factor(day, levels = sort(unique(day)))
#   )


p_composition_pre_postAbx <- ggplot(
  composition_pre_postAbx,
  aes(x = factor(day), y = relAbundance, fill = Family)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    size = 0.2
  ) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~subject) +
  labs(
    title = "Relative Abundance of Microbial Families at Passage 8",
    x = "Study Day",
    y = "Rel. Abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

savePNGPDF(paste0(OUTDIR, "compostionPrePostAbx"), p_composition_pre_postAbx, 8, 14)


### Fold-change after abx (day36/day29) by family

# Which families are the most strongly impacted by antibiotics?
  

# Filter for day 29 & 36 and only antibiotic-treated subjects
fc_pre_post_abx <- e0026 %>%
  filter(day %in% c("029", "036"), antibiotic == 1) %>%
  select(biosample1, subject, household, day, Family, OTU, relAbundance)

# Summarize relAbundance per OTU per subject per day
fc_pre_post_abx <- fc_pre_post_abx %>%
  group_by(subject, day, OTU, Family) %>%
  summarise(total_abundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = day,
    values_from = total_abundance,
    names_prefix = "relAbundance_day_"
  ) %>%
  replace_na(list(relAbundance_day_029 = 0, relAbundance_day_036 = 0)) %>%
  mutate(
    relAbundance_day_029 = if_else(relAbundance_day_029 == 0, 1e-4, relAbundance_day_029),
    relAbundance_day_036 = if_else(relAbundance_day_036 == 0, 1e-4, relAbundance_day_036),
    fold_change = relAbundance_day_036 / relAbundance_day_029,
    log2_fc = log2(fold_change)
  )


fc_pre_post_abx <- fc_pre_post_abx %>%
  filter(Family %in% top_families)

p_fc_pre_post_abx <- ggplot(fc_pre_post_abx, aes(x = fct_reorder(Family, -log2_fc), y = log2_fc, fill = Family)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_manual(values = my_colors)+
  labs(
    title = "Impact of Antibiotics on Family Relative Abundance",
    x = "Family",
    y = "Log2 Fold Change (Day 36 / Day 29)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "none"
  )

savePNGPDF(paste0(OUTDIR, "fold-changePrePostAbx"), p_fc_pre_post_abx, 6, 12)

