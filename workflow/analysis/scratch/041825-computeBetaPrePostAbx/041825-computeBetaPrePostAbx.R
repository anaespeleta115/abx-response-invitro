

### JSD between communities of day 1&29 versus JSD of communities of day 29&36

# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/041825-computeBetaPassages/041825-computeBetaPassages.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/041825-computeBetaPrePostAbx/out/"


# Filter for day combinations of interest. Maybe make a new column that has the combinations and from there decide which ones to keep.

e0026_jsd_study_days <- e0026_beta2 %>%
  mutate(day_pair = paste(day1, day2, sep = "_")) %>% 
  filter(subject1 == subject2 & day1 != day2 & passage1 == 8 & passage2 == 8, antibiotic1 == 1)
  
e0026_jsd_study_days <- e0026_jsd_study_days %>%  
  mutate(
    day_pair = map_chr(str_split(day_pair, "_"), ~ paste(sort(.x), collapse = "_"))
  )


# # Only consider day subject A samples as day 2
# e0026_beta_clean <- 
#   e0026_beta_clean %>% 
#   filter(antibiotic1 == 1)

e0026_jsd_study_days_filtered <- e0026_jsd_study_days %>% filter(day_pair %in% c("001_029", "001_036", "001_064"))

p_betaPrePostAbx <- ggplot(
  e0026_jsd_study_days_filtered,
  aes(x = day_pair, y = jsd)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.9,
    size = 1,
    fill = "#8FC0CEFF"
  # )  +
  # geom_jitter(
  #   color = "black",
  #   position = position_jitter(width = 0.1),
  #   shape = 21,
  #   stroke = 0.3,
  #   size = 1.5,
  #   alpha = 0.8
  )+
  labs(
    title = "",
    x = "Study day pair",
    y = "Jensen-Shannon Divergence"
  )  +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  DEFAULTS.THEME_PRINT



savePNGPDF(paste0(OUTDIR, "betaPrePostAbx"), p_betaPrePostAbx, 2, 2)

