source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-ploteACE010compositions/out/"


# ------------- Plot communities

# Plot compositions for day 1
eACE010_01_community_compositions <- eACE010_data %>% filter(round2plate == "A", day == "001") %>% 
  ggplot(aes(x = biosample, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~day)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT


# Plot compositions for day 29
eACE010_29_community_compositions <- eACE010_data %>% filter(round2plate == "A", day == "029") %>% 
  ggplot(aes(x = biosample, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~day)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT


# Plot compositions for day 36
eACE010_36_community_compositions <- eACE010_data %>% filter(round2plate == "A", day == "036") %>% 
  ggplot(aes(x = biosample, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~day)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT


# Plot compositions for day 64
eACE010_64_community_compositions <- eACE010_data %>% filter(round2plate == "A", day == "064") %>% 
  ggplot(aes(x = biosample, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~day)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT



final_plot <- (eACE010_01_community_compositions +
                eACE010_29_community_compositions + 
                eACE010_36_community_compositions +
                eACE010_64_community_compositions)

final_plot

savePNGPDF(paste0(OUTDIR, "eACE010_community_compositions"), final_plot, 7, 7)
