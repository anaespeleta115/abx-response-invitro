# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getRecipient-MixtureRichness/072125-getRecipient-MixtureRichness.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072225-plotRecipientRichness/out/"

# Plot
p_richness_day <- 
  ggplot(recipient_richness, aes(x = factor(day), y = species_richness)) +
  geom_line(aes(group = subject, color = subject)) +
  geom_point(aes(color = subject), size = 1.5) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "Change in Species Richness Pre- and Post-Abx",
    x = "Day",
    y = "Species Richness",
    fill = "Community Type"
  )+
  facet_wrap(~subject)+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "recipientRichness-day"), p_richness_day, 4, 8)
