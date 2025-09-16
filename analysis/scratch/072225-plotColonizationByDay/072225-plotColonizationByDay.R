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


# Statistically test 
colonization_prop_results$day <- as.factor(colonization_prop_results$day)
comparisons <- combn(levels(colonization_prop_results$day), 2, simplify = FALSE)

wilcoxon_results_colonization_prop <- left_join(
  colonization_prop_results %>%
    rstatix::wilcox_test(prop_colonizers ~ day, comparisons = comparisons),
  colonization_prop_results %>% 
    rstatix::wilcox_effsize(prop_colonizers ~ day, comparisons = comparisons),
  by = c(".y.", "group1", "group2", "n1", "n2")
) %>%
  mutate(
    comparison = paste0(group1, " vs ", group2),
    summary = paste0(
      "Wilcoxon rank-sum two-sided test, ", comparison,
      ": n=", n1 + n2,
      ", r=", round(effsize, 2),
      ", p=", signif(p, 3),
      ", adjusted p-value=", signif(p.adj, 3),
      ", ", p.adj.signif,
      ", custom adjusted p-value=", p
    )
  )


write.csv(
  wilcoxon_results_colonization_prop,
  "C:/abx-response-invitro/analysis/scratch/072225-plotColonizationByDay/out/wilcoxon_results_colonization_prop.csv",
  row.names = FALSE
)



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


# For the poster, this plot was made on replicate 2 and with a higher limit of detection (stricter cutoff)

prop_colonization_day <- 
  ggplot(colonization_prop_results, aes(x = factor(day), y = prop_colonizers, fill = recipient)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  labs(
    title = "",
    x = "Study day",
    y = "Colonization efficacy",
    fill = ""
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = PALETTE.SUBJECT)+
  facet_wrap(~recipient)+
  DEFAULTS.THEME_PRINT+
  theme(legend.position = "none")

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp-day"), prop_colonization_day, 3, 3)


prop_colonization_day_XEA <- 
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

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp-day_XEA"), prop_colonization_day_XEA, 2, 2)


