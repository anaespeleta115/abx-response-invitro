source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-plot-e0018-compositions/out/"


# mutate column to specify type of introduction (self, cohabiting, unrelated)

# plot community compositions for antibiotic-taking subjects and autologous mixtures

e0018_autologous_mixes <- e0018_data_full %>% 
  filter(mixture == "A-pre" | mixture == "A-post" | mixture == "A-pre to A-post") %>% 
  mutate(mix_type = paste(mixture, inoculationReplicate, sep = "-")) %>% 
  filter(!(mixture == "A-pre to A-post" & ratio == "1:1")) %>% 
  mutate(mix_type = fct_relevel(mix_type, "A-pre-1", "A-post-1", "A-pre to A-post-1", "A-pre to A-post-2", "A-pre to A-post-3"))


p_e0018_autologous_mixes <- e0018_autologous_mixes %>% 
  filter(household %in% c("A", "L", "P", "Q"), !mix_type %in% c("A-pre to A-post-2", "A-pre to A-post-3")) |> 
  ggplot(aes(x = mix_type, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Mixture", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~household, ncol = 11, nrow = 2)+
  DEFAULTS.THEME_PRINT

p_e0018_autologous_mixes


savePNGPDF(paste0(OUTDIR, "p_subset_e0018_autologous_mixes_non-updated"), p_e0018_autologous_mixes, 2, 2)


# plot only the autologous mixtures

p_e0018_autologous_mixes_only <- e0018_autologous_mixes %>%
  filter(!is.na(biosample2), ratio == "1:100") %>% 
  ggplot(aes(x = inoculationReplicate, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Mixture", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~household)+
  DEFAULTS.THEME_PRINT

p_e0018_autologous_mixes_only



# plot only the cohabiting mixtures


p_e0018_cohabiting_mixes <- e0018_cohabiting_mixes %>% 
  filter(mix_type == "A-pre-1" | mix_type == "A-post-1" | mix_type == "B to A-pre-1" |
           mix_type == "B to A-post-1" | mix_type == "B-1" | mix_type == "B to A-pre-2" |
           mix_type == "B to A-post-2" | mix_type == "B to A-pre-3" |
           mix_type == "B to A-post-3") %>%
  filter(!(mixture == "B to A-pre" & ratio == "1:1"), !(mixture == "B to A-post" & ratio == "1:1")) %>% 
  ggplot(aes(x = mix_type, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Mixture", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~household)+
  DEFAULTS.THEME_PRINT

p_e0018_cohabiting_mixes


savePNGPDF(paste0(OUTDIR, "p_e0018_cohabiting_mixes"), p_e0018_cohabiting_mixes, 3.5, 6)


# plot only the unrelated mixtures
# only keep one replicate of each the pre and post-mixtures

p_e0018_unrelated_mixes <- e0018_unrelated_mixes %>% 
  mutate(a_to_j = ifelse(household %in% LETTERS[1:10], "A-J", "K-V")) %>% 
  filter(a_to_j == "A-J", mix_type == "A-pre-1" | mix_type == "A-post-1" | mix_type == "C to A-pre-1" | 
                                mix_type == "C to A-post-1") %>% 
  ggplot(aes(x = mix_type, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Mixture", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  facet_wrap(~household)+
  DEFAULTS.THEME_PRINT

p_e0018_unrelated_mixes


savePNGPDF(paste0(OUTDIR, "p_e0018_unrelated_mixes"), p_e0018_unrelated_mixes, 3, 5)


