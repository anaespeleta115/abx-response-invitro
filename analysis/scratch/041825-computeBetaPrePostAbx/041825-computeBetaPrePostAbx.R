

### JSD between communities of day 1&29 versus JSD of communities of day 29&36

# Load data
source("C:/abx-response-invitro/analysis/scratch/041825-computeBetaPassages/041825-computeBetaPassages.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/041825-computeBetaPrePostAbx/out/"

# Filter for day combinations of interest. Maybe make a new column that has the combinations and from there decide which ones to keep.

e0026_beta2 <- e0026_beta2 %>%
  mutate(day_pair = paste(day1, day2, sep = "_")) %>% 
  filter(subject1 == subject2 & day1 != day2 & passage1 == 8 & passage2 == 8, antibiotic1 == 1)

# # Only consider day subject A samples as day 2
# e0026_beta_clean <- 
#   e0026_beta_clean %>% 
#   filter(antibiotic1 == 1)


p_betaPrePostAbx <- ggplot(
  e0026_beta2,
  aes(x = day_pair, y = jsd, fill = day_pair)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.9,
    size = 1
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
    title = "Divergence Across Samples Non-Abx Subjects",
    x = "Day Pair",
    y = "Jensen-Shannon Divergence"
  )  +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



savePNGPDF(paste0(OUTDIR, "betaPrePostAbx"), p_betaPrePostAbx, 6, 12)

