# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/071025-inVitro-inVivoFoldChange/out/"



### Fold-change after abx (day36/day29) by family

# Which families are the most strongly impacted by antibiotics?


# Filter for day 29 & 36 and only antibiotic-treated subjects
fc_pre_post_abx <- e0026 %>%
  filter(day %in% c("029", "036"), antibiotic == 1) %>%
  select(biosample1, subject, household, day, Family, OTU, relAbundance, passage)

# Summarize relAbundance per OTU per subject per day
fc_pre_post_abx <- fc_pre_post_abx %>%
  group_by(biosample1, day, OTU, Family, passage) %>%
  dplyr::summarise(total_abundance = sum(relAbundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = day,
    values_from = total_abundance,
    names_prefix = "relAbundance_day_"
  ) %>%
  replace_na(list(relAbundance_day_029 = 0, relAbundance_day_036 = 0)) %>%
  mutate(
    relAbundance_day_029 = if_else(relAbundance_day_029 == 0, 1e-4, relAbundance_day_029),
    relAbundance_day_036 = if_else(relAbundance_day_036 == 0, 1e-4, relAbundance_day_036),
    fold_change = relAbundance_day_036 / relAbundance_day_029,
    log10_fc = log10(fold_change)
  )

# In vitro fold change

p_fc_pre_post_abx_vitro <- fc_pre_post_abx %>% 
  filter(passage == 8) %>% 
  ggplot(aes(x = fct_reorder(Family, -log10_fc), y = log10_fc, fill = Family)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_manual(values = my_colors)+
  labs(
    title = "Impact of Antibiotics on Family Relative Abundance in Vitro (Passage 8)",
    x = "Family",
    y = "Log10 Fold Change
     (Day 36 / Day 29)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "fold-changePrePostAbx-inVitro"), p_fc_pre_post_abx_vitro, 3, 6)


# In vivo fold-change 

p_fc_pre_post_abx_vivo <- fc_pre_post_abx %>% 
  filter(passage == 0) %>% 
  ggplot(aes(x = fct_reorder(Family, -log10_fc), y = log10_fc, fill = Family)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_manual(values = my_colors)+
  labs(
    title = "Impact of Antibiotics on Family Relative Abundance in Vivo (Passage 0)",
    x = "Family",
    y = "Log10 Fold Change
     (Day 36 / Day 29)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "fold-changePrePostAbx-inVivo"), p_fc_pre_post_abx_vivo, 3, 6)


### Correlate Fold-Change

# in vitro
fc_pre_post_abx_P8 <- fc_pre_post_abx %>%
  filter(passage == 8) %>%
  group_by(Family) %>% 
  dplyr::summarise(median_fc_p8 = median(log10_fc))

# in vivo
fc_pre_post_abx_P0 <- fc_pre_post_abx %>%
  filter(passage == 0) %>%
  group_by(Family) %>% 
  dplyr::summarise(median_fc_p0 = median(log10_fc))

fc_invivo_invitro <- full_join(fc_pre_post_abx_P0, fc_pre_post_abx_P8, by = c("Family"))

# Remove rows with any missing values before calculating correlation
fc_invivo_invitro_clean <- fc_invivo_invitro %>%
  drop_na(median_fc_p0, median_fc_p8)


# Pearson correlation (assumes normal distribution)
cor_pearson <- cor(fc_invivo_invitro_clean$median_fc_p0,
                   fc_invivo_invitro_clean$median_fc_p8,
                   method = "pearson")

# Spearman correlation (rank-based, more robust to outliers)
cor_spearman <- cor(fc_invivo_invitro_clean$median_fc_p0,
                    fc_invivo_invitro_clean$median_fc_p8,
                    method = "spearman")


p_fc_invivo_invitro <- ggplot(fc_invivo_invitro_clean, aes(x = median_fc_p0, y = median_fc_p8, color = Family)) +
  geom_point(alpha = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(
    title = "",
    x = "Median Log10 Fold Change 
    Pre-Post Abx (in vivo)",
    y = "Median Log10 Fold Change 
    Pre-Post Abx (in vitro)"
  )+
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")


# Format label text
cor_label <- paste0("Pearson r = ", round(cor_pearson, 2))

# Add annotation to your existing plot
p_fc_invivo_invitro <- p_fc_invivo_invitro +
  annotate("text",
           x = min(fc_invivo_invitro_clean$median_fc_p0, na.rm = TRUE),
           y = max(fc_invivo_invitro_clean$median_fc_p8, na.rm = TRUE),
           label = cor_label,
           hjust = 0, vjust = 1,
           size = 3.5) +
  DEFAULTS.THEME_PRINT



savePNGPDF(paste0(OUTDIR, "fc_invivo_invitro"), p_fc_invivo_invitro, 3, 4)



