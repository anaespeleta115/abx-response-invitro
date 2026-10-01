source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260921-get-strain-log-info/1-260921-get-strain-log-info.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260925-get-strain-media-types/1-260925-get-strain-media-types.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260929-optimize-purity-numPlates/out/"

# In this script, I take the strain log (adjusted in the script sourced above) and filter out any low purity strains, 
# put all isolates of the same species together, and only keep the top 2 isolates in terms of their purity. After, 
# order the isolates based on the alphabetical order of their glycerol source plates.

# subset log to keep only relevant info
unique_strains <- strain_log |> 
  select(strainID, TAX_family, TAX_genus, TAX_species, TAX_purity, ISO_plate, ISO_liquid, LOC_culture, LOC_glycerol, LOC_glycerol_plateID) |> 
  group_by(TAX_species) |>
  filter(TAX_purity > 0.8) |> 
  slice_max(order_by = TAX_purity, n = 1, with_ties = TRUE) # with ties = FALSE ensures that we only keep 1 even if they have the same value.


# try ordering species by plate alphabetical order
unique_strains_ordered <- unique_strains |> 
  group_by(TAX_species) |> 
  arrange(LOC_glycerol_plateID) |> 
  slice_head(n=5)


# Add in media information:

unique_strains_final <- unique_strains_ordered |>
  left_join(species_by_media |> select(TAX_species, media_category), by = "TAX_species")


# To make final list of strains, filter out the strains with Unknown family and 
# filter out the fifth strain of common species: KSX00721, KSX00333, KSX00351, KSX00381

unique_strains_final <- unique_strains_final |> 
  filter(TAX_family != "Unknown", strainID != "KSX00721", strainID != "KSX00333", strainID != "KSX00351",
         strainID != "KSX00381", strainID != "KSX00304", strainID != "KSX00381" )

# Save this final dataset into a text file
write.table(unique_strains_final ,
            file = paste0(OUTDIR,"unique_strains_final.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

