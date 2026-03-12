source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-loade0026MetagenomicData/031126-loade0026MetagenomicData.R")



# Plot XBA compositions through time

XBA_species_composition <- metaGdata %>% 
  filter(subject == "XBA") %>% 
  ggplot(aes(x = day, y = relative_abundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))+
  DEFAULTS.THEME_PRINT

XBA_species_composition


savePNGPDF(paste0(OUTDIR, "species_abundance_XAA"), species_abundance_XAA_plot, 3, 2)