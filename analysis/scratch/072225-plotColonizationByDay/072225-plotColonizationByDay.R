# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getColonization/072125-getColonization.R")
source("C:/abx-response-invitro/analysis/scratch/072225-getColonizationProportion/072225-getColonizationProportion.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072225-plotColonizationByDay/out/"


SUBJECT <- c("XBA", "XDA", "XEA", "XKA")
PALETTE.SUBJECT <- c(  "#345995","#AD0505", "#daa520","#619E00")
names(PALETTE.SUBJECT) <- SUBJECT


total_colonization <- colonization_results %>%
  group_by(mixture, biosample1) %>%
  summarize(total_colonizers = sum(colonization)) %>%
  mutate(day = str_sub(mixture, -3), recipient = str_sub(biosample1, 0, 3))


# Plots

colonization_day <- 
  ggplot(total_colonization %>% mutate(day = as.numeric(day)-29), aes(x = factor(day), y = total_colonizers, fill = factor(day))) +
  geom_boxplot() +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Colonization proportion",
    fill = ""
  ) + 
  scale_fill_manual(values = PALETTE.DAY)+
  facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "mixtureColonization-day"), colonization_day, 4, 5)




prop_colonization_day <- 
  ggplot(colonization_prop_results %>% mutate(day = as.numeric(day)-29), aes(x = factor(day), y = prop_colonizers, fill = recipient)) +
  geom_boxplot() +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Colonization efficacy",
    fill = ""
  ) + 
  scale_fill_manual(values = PALETTE.SUBJECT)+
  facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp-day"), prop_colonization_day, 4, 4)


prop_colonization_day_XBA <- 
  ggplot(colonization_prop_results %>% filter(recipient == "XEA") %>% mutate(day = as.numeric(day)-29), aes(x = factor(day), y = prop_colonizers, fill = factor(day))) +
  geom_boxplot() +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Subject XEA
    Colonization proportion",
    fill = "Study Day"
  ) + 
  scale_fill_manual(values = PALETTE.DAY)+
  DEFAULTS.THEME_PRINT+
  theme(legend.position = "none")

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp-day_XEA"), prop_colonization_day_XBA, 2, 2)


