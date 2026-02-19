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


# Extract the top 25 families by relative abundance to make plots better
top_families <- e0026_all_passages %>%
  group_by(Family) %>%
  dplyr::summarise(total_abundance = sum(relAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Family)

# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositionsByPassage_XTA <- ggplot(e0026_all_passages %>% filter(subject == "XTA", Family %in% top_families), aes(x = factor(passage), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Passage number", y = "Relative abundance", fill = "Family") +
  theme(
    legend.position = "right"
  ) +
  facet_wrap(~ biosample1)+
  theme(legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines"), legend.text = element_text(size = 8))+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "compositionsByPassage_XTA"), p_compositionsByPassage_XTA, 1.5, 2)



legend_only <- get_legend(
  p_compositionsByPassage_XTA +
    guides(fill = guide_legend(ncol = 1)) # change 'fill' to your legend aesthetic
)

savePNGPDF(paste0(OUTDIR, "legend"), legend_only, 2.5, 1.5)
