# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getRecipient-MixtureRichness/072125-getRecipient-MixtureRichness.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072225-plotRecipientRichness/out/"

SUBJECT <- c("XBA", "XDA", "XEA", "XKA")
PALETTE.SUBJECT <- c(  "#345995","#AD0505", "#daa520","#619E00")
names(PALETTE.SUBJECT) <- SUBJECT


# Plot
p_richness_day <- 
  ggplot(recipient_richness %>% mutate(day = as.numeric(day)), aes(x = factor(day), y = species_richness, group = subject, color = subject)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(limits = c(0, NA)) + 
  labs(
    title = "",
    x = "Study day",
    y = "Species richness",
    color = "Recipient"
  ) +
  scale_color_manual(values = PALETTE.SUBJECT) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "recipientRichness-day_rep1"), p_richness_day, 1.5, 3)
