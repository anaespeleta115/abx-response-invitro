source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loadHouseholdData/031126-loadHouseholdData.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-plotHouseholdCompositions/out/"


household_data_filtered <- household_data %>% 
  filter(timepoint %in% c("29", "36", "64"), str_detect(sample, "X"))


# Plot day 1, 29, 36 and 64 compositions through time

household_species_composition <- household_data_filtered %>% 
  ggplot(aes(x = as.factor(timepoint), y = relative_abundance, fill = family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~subject)+
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT

MG_species_composition

savePNGPDF(paste0(OUTDIR, "household_species_composition"), household_species_composition, 5, 5)

