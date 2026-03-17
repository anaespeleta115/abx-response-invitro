source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-plotSpeciesRichnessMG/031726-plotSpeciesRichnessMG.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-plotSpeciesRichnessHousehold/031726-plotSpeciesRichnessHousehold.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-compareRichnessp0p8/out/"




PASSAGE <- c(0, 8)
PALETTE.PASSAGE <- c("gray80","#7570B3")
names(PALETTE.PASSAGE) <- PASSAGE



#join together the two species richness datasets (household and in vitro)

MG_p8_species_richness <- MG_p8_species_richness %>% 
  mutate(passage = "8")

MG_p0_species_richness <- MG_p0_species_richness %>% 
  mutate(passage = "0")

invivo_invitro_richness <- rbind(MG_p0_species_richness, MG_p8_species_richness) 
  


# plot only the antibiotic taking subjects, but now color by passage
p_invivo_invitro_richness <- invivo_invitro_richness %>% mutate(day = as.numeric(day)-29) %>% 
  ggplot(aes(x = factor(day), y = species_richness, fill = factor(passage))
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
    fill = "Passage"
  ) +
  ylim(0, 80)+
  scale_fill_manual(values = PALETTE.PASSAGE)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 0.5, size = 8),
    axis.text.y = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8)
  )

savePNGPDF(paste0(OUTDIR, "invivo_invitro_richness"), p_invivo_invitro_richness, 2, 3)
