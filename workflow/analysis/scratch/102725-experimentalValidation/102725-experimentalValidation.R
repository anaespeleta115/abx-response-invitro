source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102725-experimentalValidation/out/"


# --------------------------- PLOT XKB COMPOSITION Passsage 0 & passage 8 -------------------------

p_XKB_compositions <- ggplot(e0026 %>% filter(biosample1 == "XKB-029", passage %in% c(0, 8)), aes(x = biosample1, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  facet_wrap(~passage)+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "KXB_composition"), p_XKB_compositions, 3, 3)
