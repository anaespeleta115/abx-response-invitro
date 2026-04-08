source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-ploteAME004compositions/out/"



# Plot mixture communities
eAME004_mixture_compositions <- eAME004_data %>% filter(replicate == 1, subject == "XGA") %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "mixture", y = "relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eAME004_mixture_compositions


savePNGPDF(paste0(OUTDIR, "eAME004_community_compositions"), eAME004_mixture_compositions, 3, 5)


eAME004_recipient_compositions <- eAME004_data %>% filter(replicate == 1, biosample2 == "blank") %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "mixture", y = "relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eAME004_recipient_compositions