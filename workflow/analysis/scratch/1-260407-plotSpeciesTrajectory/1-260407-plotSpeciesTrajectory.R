source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-ploteACE010compositions/1-260406-ploteACE010compositions.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-geteAME004Richness/1-260408-geteAME004SpeciesRichness.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260407-plotSpeciesTrajectory/out/"


subjectsMixed <- c("XEA", "XDA", "XCA", "XFA", "XGA", "XSA", "XJA", "XLA")

# Plot eACE010 species richness through time
eACE010_richness_by_day <- richness %>% filter(subject %in% subjectsMixed) %>%  
  ggplot(aes(x = factor(day), y = nASVs))+
  geom_line(aes(group=1)) +
  geom_point(size = 0.5)+
  # geom_text(
  #   data = richness,
  #   aes(x = day, y = y_pos, label = nASVs),
  #   inherit.aes = FALSE,
  #   size = 2
  # ) +
  labs(x = "study day", y = "species richness")+
  facet_wrap(~subject)+
  ylim(0, 70)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eACE010_richness_by_day

savePNGPDF(paste0(OUTDIR, "eACE010_richness_by_day_subjectsMixed"), eACE010_richness_by_day, 2.5, 2.5)

### PLOT EAME004 SPECIES RICHNESS IN THE SAME MANNER, COMPARE THE TWO FIGURES.

eAME004_richness <- eAME004_recipient_diversity %>% 
  ggplot(aes(x = factor(day), y = nASVs))+
  geom_line(aes(group=1)) +
  geom_point(size = 0.5)+
  # geom_text(
  #   data = eAME004_delta_diversity,
  #   aes(x = day, y = y_pos, label = nASVs),
  #   inherit.aes = FALSE,
  #   size = 2
  # ) +
  labs(x = "study day", y = "species richness")+
  facet_wrap(~subject)+
  ylim(0, 70)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eAME004_richness

savePNGPDF(paste0(OUTDIR, "eAME004_richness"), eAME004_richness, 2.5, 2.5)
