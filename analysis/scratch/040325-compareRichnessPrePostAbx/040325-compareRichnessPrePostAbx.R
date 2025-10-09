# Species Richness: calculate number of OTU's with a relative abundance > 0.1%. 
# Plot the distribution of species richness at each timepoint, for abx- and non-abx subjects.

# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/040325-comparerichnessPrePostAbx/out/"

# Set palette
pal <- pnw_palette("Sunset2", 2, type = "discrete")

ABX <- c("No", "Yes")
PALETTE.ABX <- c("gray80","#88CCEE")
names(PALETTE.ABX) <- ABX


PASSAGE <- c(0, 8)
PALETTE.PASSAGE <- c(  "gray80", "#B084CC")
names(PALETTE.PASSAGE) <- PASSAGE


e0026_richness <- e0026_richness %>%
  ungroup() %>%  
  mutate(
    day = as.character(day),
    antibiotic = factor(antibiotic, levels = c(0, 1), labels = c("No", "Yes")),
    abxDay = factor(paste0(antibiotic, "_", day)),
    species_richness = as.numeric(species_richness)
  ) %>%
  filter(passage %in% c(0, 8)) %>%
  droplevels()

# Statistically test 

comparisons <- combn(levels(e0026_richness$abxDay), 2, simplify = FALSE)
# 
# wilcoxon_results_p0 <- left_join(
#   e0026_richness %>%
#     filter(passage == 0) %>% 
#     rstatix::wilcox_test(species_richness ~ abxDay, comparisons = comparisons),
#   e0026_richness %>%
#     filter(passage == 0) %>% 
#     rstatix::wilcox_effsize(species_richness ~ abxDay, comparisons = comparisons),
#   by = c(".y.", "group1", "group2", "n1", "n2")
# ) %>%
#   mutate(
#     comparison = paste0(group1, " vs ", group2),
#     summary = paste0(
#       "Wilcoxon rank-sum two-sided test, ", comparison,
#       ": n=", n1 + n2,
#       ", r=", round(effsize, 2),
#       ", p=", signif(p, 3),
#       ", adjusted p-value=", signif(p.adj, 3),
#       ", ", p.adj.signif,
#       ", custom adjusted p-value=", p
#     )
#   )


wilcoxon_results_p8 <- left_join(
  e0026_richness %>%
    filter(passage == 8) %>% 
    rstatix::wilcox_test(species_richness ~ abxDay, comparisons = comparisons),
  e0026_richness %>%
    filter(passage == 8) %>% 
    rstatix::wilcox_effsize(species_richness ~ abxDay, comparisons = comparisons),
  by = c(".y.", "group1", "group2", "n1", "n2")
) %>%
  mutate(
    comparison = paste0(group1, " vs ", group2),
    summary = paste0(
      "Wilcoxon rank-sum two-sided test, ", comparison,
      ": n=", n1 + n2,
      ", r=", round(effsize, 2),
      ", p=", signif(p, 3),
      ", adjusted p-value=", signif(p.adj, 3),
      ", ", p.adj.signif,
      ", custom adjusted p-value=", 2 * p
    )
  )

valid_comparisons <- c("Yes_029 vs Yes_036", "Yes_029 vs Yes_064")

# wilcoxon_results_p0 <- wilcoxon_results_p0 %>%
#   filter(comparison %in% valid_comparisons)
wilcoxon_results_p8 <- wilcoxon_results_p8 %>%
  filter(comparison %in% valid_comparisons)


write.csv(
  wilcoxon_results_p8,
  "/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/040325-comparerichnessPrePostAbx/out/wilcoxon_results_richness.csv",
  row.names = FALSE
)



e0026_richness_filtered <- e0026_richness %>% 
  filter(antibiotic != "0", passage %in% c(0,8)) %>% 
  mutate(
    passage = factor(passage, levels = c(0, 8))
  )


p_richness_time_filtered <- ggplot(
  e0026_richness_filtered %>% mutate(day = as.numeric(day)-29),
  aes(x = factor(day), y = species_richness, fill = passage)
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
    title = "",
    x = "Study day",
    y = "Species richness",
    fill = "Passage"
  ) +
  scale_fill_manual(values=PALETTE.PASSAGE)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 0.5, size = 8),
    axis.text.y = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8)
  )

savePNGPDF(paste0(OUTDIR, "richnessByTime_filtered"), p_richness_time_filtered, 2, 3)
