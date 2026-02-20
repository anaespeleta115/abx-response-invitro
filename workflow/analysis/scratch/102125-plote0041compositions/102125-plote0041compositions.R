# load scripts
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102625-gete0041LostTaxa/102625gete0041LostTaxa.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-plote0041compositions/out/"


# --------------------------- PLOT DONOR/RECIPIENTS -------------------------

p_donor_compositions <- ggplot(e0041_control_donors, aes(x = recipient, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  facet_wrap(~donor)+
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_donors"), p_donor_compositions, 4, 4)

p_recipient_compositions <- ggplot(e0041_control_recipients %>% filter(recipient == "post-abx-V2"), aes(x = recipient, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  facet_wrap(~donor)+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_recipients_V2"), p_recipient_compositions, 2.25, 1)

# --------------------------- PLOT PRE-ABX MIXTURES ---------------------------------

p_pre_abx <- ggplot(mixture_colonization_full %>% filter(recipient == "pre-abx", donor != "XBB-029" & donor != "XEB-029" & donor != "XDB-029"), aes(x = recipient, y = relAbundance, fill = Family, alpha = factor(colonizer))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  facet_wrap(~donor)+
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_pre-abx_colonization"), p_pre_abx, 4, 4)

# --------------------------- PLOT PRE-ABX MIXTURES (with different labels for colonization) ---------------------------------


p_post_abx_v1 <- ggplot(mixture_colonization_full %>% filter(recipient == "post-abx-V1", donor != "XBB-029" & donor != "XEB-029" & donor != "XDB-029"), aes(x = recipient, y = relAbundance, fill = Family, alpha = factor(diff_colonizer_v1))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  facet_wrap(~donor)+
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_post-abx-V1_DIFFcolonization"), p_post_abx_v1, 4, 4)

p_post_abx_v2 <- ggplot(mixture_colonization_full %>% filter(recipient == "post-abx-V2", donor == "XIB-029"), aes(x = recipient, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  # facet_wrap(~donor)+
  # theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_post-abx-V2_DIFFcolonization"), p_post_abx_v2, 2, 1)


# --------------------------------- PLOT ALL COMPOSITIONS BY WELL ---------------------------------------
plot_list <- list()

for (i in seq_along(wells)) {
  w <- wells[i]
  
  d <- e0041.A %>% 
    filter(.data$well == w)
  
  p <- ggplot(d, aes(x = 1, y = relAbundance, fill = Family)) +
    geom_col(color = "black", linewidth = 0.1) +
    scale_fill_manual(values = my_colors) +
    labs(title = w, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    DEFAULTS.THEME_PRINT
  
  plot_list[[i]] <- p
}

# Adjust ncol to plate/grid layout (12 for a 96-well plate row)
combined_plot <- wrap_plots(plot_list, ncol = 12)
combined_plot

savePNGPDF(paste0(OUTDIR, "all_compositions"), combined_plot, 10, 10)

# --------------------------------- PLOT PRE-ABX LOST TAXA V1 ---------------------------------------


p_pre_abx_lost_v1 <- ggplot(e0041_control_recipients %>% filter(recipient == "pre-abx") %>% mutate(lost_V1 = ifelse(OTU %in% post_abx_lost_taxa_v1, 1, 0)), aes(x = recipient, y = relAbundance, fill = Family, alpha = factor(lost_V1))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_pre-abx_lost_V1"), p_pre_abx_lost_v1, 2, 1)




# --------------------------------- PLOT PRE-ABX LOST TAXA V2 ---------------------------------------

# , alpha = factor(lost_V2))

p_pre_abx_lost_v2 <- ggplot(e0041_control_recipients %>% filter(recipient == "pre-abx") %>% mutate(lost_V2 = ifelse(OTU %in% post_abx_lost_taxa_v2, 1, 0)), aes(x = recipient, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_pre-abx_lost_V2"), p_pre_abx_lost_v2, 2, 1)


# Q: ARE ALL LOST TAXA GAINED BACK IN THE MIXTURE? SUTTERELLACEAE AND ENTEROBACTERIACEAE ARE GAINED BACK IN ALL MIXTURES.


