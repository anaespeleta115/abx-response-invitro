### Compositional changes across passages

# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/040325-plotCompositions/out/"

# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositionsByPassage <- ggplot(e0026_all_passages, aes(x = factor(passage), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Passage number", y = "Relative abundance", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  facet_wrap(~ biosample1)+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "compositionsByPassage"), p_compositionsByPassage, 3, 4)



# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositionsByPassage_XTA <- ggplot(e0026_all_passages %>% filter(subject == "XTA"), aes(x = factor(passage), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Passage number", y = "Relative abundance", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  facet_wrap(~ biosample1)+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "compositionsByPassage_XTA"), p_compositionsByPassage_XTA, 1.5, 2)