library(cowplot)
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loade0026MetagenomicData/031126-loade0026MetagenomicData.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-plotMGcompositions/out/"


# Plot XBA compositions through time

MG_species_composition <- metaGdata %>% 
  ggplot(aes(x = day, y = relative_abundance, fill = family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~subject)+
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

MG_species_composition

savePNGPDF(paste0(OUTDIR, "MG_species_composition"), MG_species_composition, 5, 5)


# get legend
# 
# p_with_legend <- metaGdata %>% 
#   ggplot(aes(x = day, y = relative_abundance, fill = family)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = my_colors) +
#   facet_wrap(~subject) +
#   labs(x = "Study day", y = "Relative abundance") +
#   DEFAULTS.THEME_PRINT
# 
# 
# 
# legend <- cowplot::get_legend(
#   p_with_legend + theme(legend.position = "right")
# )
# 
# savePNGPDF(paste0(OUTDIR, "metagenomic_legend"), legend, 3, 3)
