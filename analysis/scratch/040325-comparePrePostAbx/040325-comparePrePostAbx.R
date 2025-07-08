# Species Richness: calculate number of OTU's with a relative abundance > 0.1%. 
# Plot the distribution of species richness at each timepoint, for abx- and non-abx subjects.

# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/040325-comparePrePostAbx/out/"

# Combine all the day datasets

e0026_richness <- e0026_richness %>%
  mutate(
    day = factor(day, levels = sort(unique(day))),
    antibiotic = factor(antibiotic, levels = c(0, 1), labels = c("No", "Yes"))
  )


p_richness_time <- ggplot(
  e0026_richness,
  aes(x = day, y = species_richness, fill = antibiotic)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.9,
    size = 1
  )  +
  # geom_jitter(
  #   color = "black",
  #   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
  #   shape = 21,
  #   stroke = 0.3,
  #   size = 1.5,
  #   alpha = 0.8
  # ) +
  labs(
    title = "Change in Species Richness Pre- and Post Abx",
    x = "Day",
    y = "Species Richness",
    fill = "Antibiotic"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

savePNGPDF(paste0(OUTDIR, "richnessByTime"), p_richness_time, 6, 12)
