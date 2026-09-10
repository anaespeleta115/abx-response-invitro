source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260824-get-e0018-species-richness/out/"


# ---------- get species richness for recipients and donors and AUTOLOGOUS mixes


autologous_richness <- e0018_autologous_mixes %>% 
  filter(relAbundance > limit_of_detection) %>% 
  group_by(household, mix_type, biosample_full) %>% 
  summarize(richness = n())

# plot richness by 
p_autologous_richness <- autologous_richness %>% 
  ggplot(aes(x = mix_type, y = richness, color = mix_type))+
    geom_point()+
    labs(
      x = "Community",
      y = "Species Richness"
    )+
  ylim(0, 70)+
    DEFAULTS.THEME_PRINT+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )+
  facet_wrap(~ household)+
  scale_color_paletteer_d("rcartocolor::BurgYl")

p_autologous_richness

savePNGPDF(paste0(OUTDIR, "p_autologous_richness"), p_autologous_richness, 4.5, 4.5)


# ----------- get species richness for COHABITING mixtures


cohabiting_richness <- e0018_cohabiting_mixes %>% 
  filter(relAbundance > limit_of_detection) %>% 
  filter(!(mixture == "B to A-pre" & ratio == "1:1"), !(mixture == "B to A-post" & ratio == "1:1")) %>% 
  group_by(household, mix_type, biosample_full) %>% 
  summarize(richness = n())

# plot richness
p_cohabiting_richness <- cohabiting_richness %>% 
  ggplot(aes(x = mix_type, y = richness, color = mix_type))+
  geom_point()+
  labs(
    x = "Community",
    y = "Species Richness"
  )+
  ylim(0, 80)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~ household)+
  scale_color_paletteer_d("ggthemes::Classic_10_Light")

p_cohabiting_richness

savePNGPDF(paste0(OUTDIR, "p_cohabiting_richness"), p_cohabiting_richness, 4.5, 5)


# ----------- get species richness for UNRELATED mixtures



unrelated_richness <- e0018_unrelated_mixes %>% 
  mutate(a_to_j = ifelse(household %in% LETTERS[1:10], "A-J", "K-V")) %>% 
  filter(relAbundance > limit_of_detection, a_to_j == "A-J") %>%
  group_by(household, mix_type, biosample_full) %>% 
  summarize(richness = n())

# plot richness 
p_unrelated_richness <- unrelated_richness %>% 
  ggplot(aes(x = mix_type, y = richness, color = mix_type))+
  geom_point()+
  labs(
    x = "Community",
    y = "Species Richness"
  )+
  ylim(0, 80)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~ household)+
  scale_color_paletteer_d("MetBrewer::Hiroshige")

p_unrelated_richness

savePNGPDF(paste0(OUTDIR, "p_unrelated_richness"), p_unrelated_richness, 3.5, 4)

