source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031126-plotMGcompositions/031126-plotMGcompositions.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031226-plotHouseholdCompositions/031226-plotHouseholdCompositions.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/031726-compareXBAcompositions_MG/out/"


# select relevant columns for both datasets
household_filtered <- household_filtered %>% 
  select(sample, species_id, relative_abundance, subject, kingdom, phylum, class, order, family, genus, species,
         antibiotic, day) %>% 
  mutate(passage = 0)

metaGdata_filtered <- metaGdata_filtered %>% 
  select(sample, species_id, relative_abundance, subject, kingdom, phylum, class, order, family, genus, species,
         antibiotic, day) %>% 
  mutate(passage = 8)


# bind together household and invitro data
invivo_invitro_data <- rbind(household_filtered, metaGdata_filtered)

# Plot in vivo and in vitro compositions for XBA
p_compositions_invivo_invitro_XEA <- invivo_invitro_data %>% 
  filter(subject == "XKA") %>% mutate(day = as.numeric(day)-29) %>% 
  ggplot(
    aes(x = factor(day), y = relative_abundance, fill = family)
  ) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    size = 0.2
  ) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~passage) +
  labs(
    title = "",
    x = "Study day",
    y = "Relative abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none"
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XKA_MG_PrePostCompositions"), p_compositions_invivo_invitro_XEA, 3, 3)



