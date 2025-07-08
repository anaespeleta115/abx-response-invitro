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
  summarise(
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
  mutate(log2_fc = log2(fold_change))
  

p_P0P8 <- P0P8 %>%
  ggplot(aes(x = fct_reorder(Family, -log2_fc), y = log2_fc, fill = Family)) +
  geom_boxplot(width = 0.9, outlier.shape = NA) +
  scale_fill_manual(values = my_colors)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(
    title = "Distribution of Log2 Fold Change per Subject (P8 vs P0)",
    x = "Subject (biosample1)",
    y = "Log2 Fold Change"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")


savePNGPDF(paste0(OUTDIR, "compareP0P8"), p_P0P8, 6, 12)





