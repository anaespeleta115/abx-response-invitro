# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getColonization/072125-getColonization.R")
source("C:/abx-response-invitro/analysis/scratch/072225-getColonizationProportion/072225-getColonizationProportion.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072525-plotColonizationHeatmap/out/"


p_colonization_heatmap <- ggplot(colonization_prop_results, aes(x = biosample2, y = biosample1, fill = prop_colonizers)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A") +
  labs(x = "Donor", y = "Recipient", fill = "Colonization efficacy",
       title = "")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  DEFAULTS.THEME_PRINT



savePNGPDF(paste0(OUTDIR, "propColonizationHeatmap"), p_colonization_heatmap, 2, 2)

