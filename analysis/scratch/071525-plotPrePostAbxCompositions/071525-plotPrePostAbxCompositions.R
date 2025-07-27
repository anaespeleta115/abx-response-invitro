# Plot stacked bar plots with a facet for each subject, showing compositional changes over time across timepoints, at passage 8.

# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/071525-plotPrePostAbxCompositions/out/"


### Compositional changes pre and post-abx

P8_compositions <- combined_day_data %>% 
  filter(passage == 8) %>% 
  filter(subject %in% lastingResponses)

# Summarize relative abundance per biosample1 by Family
composition_pre_postAbx <- P8_compositions %>%
  group_by(biosample1, Family) %>%
  dplyr::summarise(relAbundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
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
    title = "In vitro compositions",
    x = "Study day",
    y = "Relative abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compostionPrePostAbx-P8"), p_composition_pre_postAbx, 4, 8)


### Compositional changes pre and post-abx

P0_compositions <- combined_day_data %>% 
  filter(passage == 0) %>% 
  filter(subject %in% lastingResponses)

# Summarize relative abundance per biosample1 by Family
composition_pre_postAbx0 <- P0_compositions %>%
  group_by(biosample1, Family) %>%
  dplyr::summarise(relAbundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    P8_compositions %>% distinct(biosample1, day, subject, household, antibiotic),
    by = "biosample1"
  ) 
# 
# %>%
#   mutate(
#     day = factor(day, levels = sort(unique(day)))
#   )


p_composition_pre_postAbx0 <- ggplot(
  composition_pre_postAbx0,
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
    title = "In vivo compositions",
    x = "Study day",
    y = "Relative abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compostionPrePostAbx-P0"), p_composition_pre_postAbx0, 4, 8)



XBA_compositions <- combined_day_data %>% 
  filter(subject == "XBA")

# Summarize relative abundance per biosample1 by Family
composition_pre_postAbx_XBA <- XBA_compositions %>%
  group_by(passage, biosample1, Family) %>%
  dplyr::summarise(relAbundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    combined_day_data %>% distinct(biosample1, day, subject, household, antibiotic) %>% filter(subject == "XBA"),
    by = "biosample1"
  ) 
# 
# %>%
#   mutate(
#     day = factor(day, levels = sort(unique(day)))
#   )


p_composition_pre_postAbx_XBA <- ggplot(
  composition_pre_postAbx_XBA,
  aes(x = factor(day), y = relAbundance, fill = Family)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    size = 0.2
  ) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~passage) +
  labs(
    title = "",
    x = "Study day",
    y = "Relative abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "compostionPrePostAbx-XBA"), p_composition_pre_postAbx_XBA, 2, 3.5)



