source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loade0026MetagenomicData/031126-loade0026MetagenomicData.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-plotSpeciesRichnessMG/out/"



ABX <- c("0", "1")
PALETTE.ABX <- c("gray80","#88CCEE")
names(PALETTE.ABX) <- ABX



MG_p8_species_richness <- metaGdata_filtered %>% 
  group_by(subject, day, antibiotic) %>% 
  summarize(species_richness = n())



MG_p8_species_richness_plot <- MG_p8_species_richness %>% mutate(day = as.numeric(day)-29) %>% 
  ggplot(aes(x = factor(day), y = species_richness, fill = factor(antibiotic))
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.9,
    size = 1
  )  +
  labs(
    title = "",
    x = "Study day",
    y = "Species richness",
    fill = "Antibiotic"
  )+
  scale_fill_manual(values = PALETTE.ABX)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 0.5, size = 8),
    axis.text.y = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8)
  )

savePNGPDF(paste0(OUTDIR, "MG_richness_time_p8"), MG_p8_species_richness_plot, 3, 4)

