source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loadHouseholdData/031126-loadHouseholdData.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-plotHouseholdCompositions/out/"



# Plot day 2, 29, 36 and 64 compositions through time

household_species_composition <- household_filtered %>% 
  ggplot(aes(x = as.factor(day), y = relative_abundance, fill = family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~subject)+
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT



savePNGPDF(paste0(OUTDIR, "household_species_composition"), household_species_composition, 5, 5)

