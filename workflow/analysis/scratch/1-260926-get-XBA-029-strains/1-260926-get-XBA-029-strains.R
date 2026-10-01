source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260921-get-strain-log-info/1-260921-get-strain-log-info.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260926-get-XBA-029-strains/out/"


XBA_029_new <- read_csv("~/Documents/GitHub/abx-response-invitro/data/C5GM3F_1_eAME006_1_taxonomy-table.csv") |> 
  mutate(subject = "XBA-029") |> 
  mutate(Family = ifelse(Genus == "Flintibacter", "Oscillospiraceae", Family))
         # Family = ifelse(Genus == "Fenollaria", "Peptoniphilaceae", Family))


# plot compositions of new XBA community
p_XBA_compositions <- XBA_029_new |> 
  ggplot(aes(x = subject, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT
  
p_XBA_compositions

savePNGPDF(paste0(OUTDIR, "p_XBA_compositions"), p_XBA_compositions, 3, 3)



# join in the strains from strain log and XBA to find exact matches.

unique_strains_joined <- unique_strains |> 
  full_join(XBA_029_new, by = c("TAX_species" = "Species", "TAX_family" = "Family"))


# Save dataset into a text file
write.table(unique_strains_joined ,
            file = paste0(OUTDIR,"unique_strains_joined.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

