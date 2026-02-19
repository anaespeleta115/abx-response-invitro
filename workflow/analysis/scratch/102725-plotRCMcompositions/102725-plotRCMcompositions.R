source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102025-loade0040data/102025-loade0040data.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102725-plotRCMcompositions/out/"



# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositions_RCM <- ggplot(e0040_RCM, aes(x = well, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color= "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Pre-abx RCM " , x = "", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT


# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositions_mBHI <- ggplot(e0040 %>% filter(community == "XEA-pre-abx"), aes(x = well, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Pre-abx mBHI" , x = "", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

combined_medias <- p_compositions_RCM / p_compositions_mBHI

savePNGPDF(paste0(OUTDIR, "combined_replicates"), combined_medias, 5, 5)
