source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-ploteAME004colonization/out/"


# set palette
DAY <- c("0", "7", "35")
PALETTE.DAY <- c(paletteer_d("nationalparkcolors::Voyageurs"),1)
names(PALETTE.DAY) <- DAY

# ------------- PLOT RAW COLONIZER COUNTS BY DAY -----------------

total_colonization <- colonization_results %>%
  group_by(mixture, biosample1) %>%
  summarize(total_colonizers = sum(colonizer)) %>% 
  mutate(
    subject1 = str_sub(biosample1, 1, -5),
    subject2 = str_sub(mixture, 1, 3),
    day = str_sub(biosample1, -3))

colonization_day <- 
  ggplot(total_colonization %>%  mutate(day = (as.numeric(as.character(day)))-29), aes(x = factor(day), y = total_colonizers)) +
  geom_boxplot(fill = "#5495CFFF", outlier.size = 0.5, ) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Number of colonizers",
    fill = ""
  ) + 
  # scale_fill_manual(values = "#5495")+
  # facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "mixtureColonization-day"), colonization_day, 2, 2)


# ------------ PLOT COLONIZATION PROPORTION BY DAY -----------------

prop_colonization_day <- 
  ggplot(colonization_prop_results %>% mutate(day = (as.numeric(as.character(day)))-29), aes(x = factor(day), y = prop_colonizers)) +
  geom_boxplot(outlier.size = 0.5, fill = "#5495CFFF") +
  # geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Colonization efficacy",
    fill = ""
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  # facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT+
  theme(legend.position = "none")

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp-day"), prop_colonization_day, 3, 4)
