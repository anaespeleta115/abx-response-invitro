source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260921-get-strain-log-info/1-260921-get-strain-log-info.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260925-get-strain-media-types/out/"


## In this script I attempt to find out which strains are capable of growing in mBHI and which can only grow on RCM.

## Brute-force way to do this:

# divide dataset into the two media types (mBHI and RCM)
mBHI <- strain_log |> 
  filter(ISO_liquid == "mBHI")
# specify media type in column names
colnames(mBHI) <- paste0("mBHI_", colnames(mBHI))


RCM <- strain_log |> 
  filter(ISO_liquid == "RCM")
# specify media type in column names
colnames(RCM) <- paste0("RCM_", colnames(RCM))

# the strains that don't match any of the mBHI strains are the ones that can only grow on RCM
media_joined <- mBHI |> select(mBHI_strainID, mBHI_TAX_species, mBHI_TAX_family, mBHI_ISO_liquid) |> 
  full_join(RCM |> select(RCM_strainID, RCM_TAX_species, RCM_TAX_family, RCM_ISO_liquid), by = c("mBHI_TAX_species" = "RCM_TAX_species"),
            relationship = "many-to-many")


# A more effective way to do this:
species_by_media <- strain_log |>
  group_by(TAX_species) |>
  summarise(
    in_mBHI = any(ISO_liquid == "mBHI"),
    in_RCM  = any(ISO_liquid == "RCM"),
    .groups = "drop"
  ) |>
  mutate(
    media_category = case_when(
      in_mBHI & in_RCM  ~ "Both",
      in_mBHI            ~ "mBHI only",
      in_RCM             ~ "RCM only"
    )
  )


# # Save dataset into a text file
# write.table(species_by_media ,
#             file = paste0(OUTDIR,"species_by_media.txt"),
#             sep = "\t",
#             row.names = FALSE,
#             quote = FALSE)



# add media type classification to unique strains

unique_strains <- unique_strains |>
  left_join(species_by_media |> select(TAX_species, media_category), by = "TAX_species")

# # Save dataset into a text file
# write.table(unique_strains ,
#             file = paste0(OUTDIR,"unique_strains.txt"),
#             sep = "\t",
#             row.names = FALSE,
#             quote = FALSE)


