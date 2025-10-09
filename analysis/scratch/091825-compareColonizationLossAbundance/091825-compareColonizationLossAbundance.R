# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/091825-compareColonizationLossAbundance/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_loss <- recipient_ASVs %>%
  filter(day == "029" | day == "036", replicate == curr_replicate, relAbundance > 1e-4) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(day, Family, subject1) %>%
  summarise(fam_relAbundance = sum(relAbundance))

recipient_fam_loss <- recipient_fam_loss %>%
  pivot_wider(
    names_from = day,
    values_from = fam_relAbundance,
    names_prefix = "recipient_fam_relAbundance_day_"
  )%>%
  replace_na(list(recipient_fam_relAbundance_day_029 = 0, recipient_fam_relAbundance_day_036 = 0)) %>%
  filter(recipient_fam_relAbundance_day_029 != 0 & recipient_fam_relAbundance_day_036 != 0) %>%  # filter to make sure that we only have taxa present in that community
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(recipient_fam_relAbundance_day_029 == 0, 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(recipient_fam_relAbundance_day_036 == 0, 1e-4, recipient_fam_relAbundance_day_036),
    recipient_ratio = recipient_fam_relAbundance_day_029 / recipient_fam_relAbundance_day_036,
    recipient_log10_ratio = log10(recipient_ratio)
  )
