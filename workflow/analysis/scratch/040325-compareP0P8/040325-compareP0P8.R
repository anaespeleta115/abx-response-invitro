# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/040325-compareP0P8/out/"


### Fold-change of each family in passage 0 versus passage 8 (p0/p8 for each subject in day 29)


# Get passage 0 and passage 8 communities for day 29

P0P8 <- e0026_day29 %>% 
  filter(passage %in% c(0, 8)) %>% 
  select(biosample1, OTU, count, passage, relAbundance, Family) %>% 
  group_by(passage, biosample1, Family) %>% 
  dplyr::summarise(
    total_abundance = sum(relAbundance, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = passage,
    values_from = total_abundance,
    names_prefix = "relAbundance_passage_"
  )%>%
  replace_na(list(relAbundance_passage_0 = 0, relAbundance_passage_8 = 0)) %>% 
  mutate(
    relAbundance_passage_0 = if_else(relAbundance_passage_0 == 0, 1e-4, relAbundance_passage_0),
    relAbundance_passage_8 = if_else(relAbundance_passage_8 == 0, 1e-4, relAbundance_passage_8)
  ) %>% 
  mutate(fold_change = relAbundance_passage_8/relAbundance_passage_0) %>% 
  filter(Family %in% top_families) %>% 
  mutate(log10_fc = log10(fold_change))
  

p_P0P8 <- P0P8 %>%
  ggplot(aes(x = fct_reorder(Family, -log10_fc), y = log10_fc, fill = Family)) +
  geom_boxplot(width = 0.9, outlier.shape = NA) +
  scale_fill_manual(values = my_colors)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred")+ 
  labs(
    title = "",
    x = "Family",
    y = "Log10 Fold Change in vitro/in vivo"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "compareP0P8"), p_P0P8, 4, 6)





