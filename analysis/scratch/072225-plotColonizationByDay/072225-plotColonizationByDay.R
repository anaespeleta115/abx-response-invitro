# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getColonization/072125-getColonization.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072225-plotColonizationByDay/out/"



total_colonization <- colonization_results %>% 
  group_by(mixture, biosample1) %>% 
  summarize(total_colonizers = sum(colonization)) %>% 
  mutate(day = str_sub(mixture, -2), recipient = str_sub(biosample1, 0, 3))


# Plot
p_colonization_day <- 
  ggplot(total_colonization, aes(x = factor(day), y = total_colonizers, fill = factor(day))) +
  geom_boxplot() +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Number of colonizers",
    fill = ""
  ) + 
  facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "mixtureColonization-day"), p_colonization_day, 3, 4)
