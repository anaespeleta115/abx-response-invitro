# load scripts
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/111425-loade0058data/111425-loade0058data.R")

source("~/Documents/GitHub/abx-response-invitro/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/111623-plotE0058compositions/out/"



# Plot the two full communities
p_full_community <- ggplot(e0058 %>% filter(community == "full_community" | community == "full_community_ss"), aes(x = community, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "p_full_community"), combined_plot, 4, 4)


# Plot the communities with bacteroides dropped
p_dropped_bacter <- ggplot(e0058 %>% filter(str_detect(community, "bacteroides")), aes(x = community, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


combined_plot <- p_full_community / p_dropped_bacter 

savePNGPDF(paste0(OUTDIR, "all_communities"), combined_plot, 4, 4)


# ------------------------------------- PLOT ALL COMMUNITIES -----------------------------------------
plot_list <- list()

for (i in seq_along(wells)) {
  w <- wells[i]
  
  d <- e0058 %>% 
    filter(.data$well == w)
  
  p <- ggplot(d, aes(x = 1, y = relAbundance, fill = Family)) +
    geom_col(color = "black", linewidth = 0.1) +
    scale_fill_manual(values = my_colors) +
    labs(title = w, x = NULL, y = NULL) +
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

savePNGPDF(paste0(OUTDIR, "all_compositions"), combined_plot, 8, 8)


