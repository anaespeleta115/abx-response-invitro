source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-ploteAME004compositions/out/"


# Plot mixture communities
eAME004_all_compositions <- eAME004_data %>% filter(replicate == curr_replicate, subject == "XJA") %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "recipient", y = "relative abundance") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  # facet_wrap(~biosample2)+
  DEFAULTS.THEME_PRINT

eAME004_all_compositions

savePNGPDF(paste0(OUTDIR, "eAME004_XJA_compositions"), eAME004_all_compositions, 5, 6)




eAME004_recipient_compositions <- eAME004_data %>% filter(replicate == curr_replicate, biosample2 == "blank") %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(x = "recipient", y = "relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eAME004_recipient_compositions

savePNGPDF(paste0(OUTDIR, "eAME004_recipient_compositions"), eAME004_recipient_compositions, 3, 4)

eAME004_donor_compositions <- eACE010_data_unpooled %>% filter(biosample %in% donor_communities) %>% 
  ggplot(aes(x = biosample, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "donor", y = "relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eAME004_donor_compositions

savePNGPDF(paste0(OUTDIR, "eAME004_donor_compositions"), eAME004_donor_compositions, 3, 4)


# put both recipient and donor compositions together
eAME004_joint_compositions <- eAME004_recipient_compositions / eAME004_donor_compositions

savePNGPDF(paste0(OUTDIR, "eAME004_joint_compositions"), eAME004_joint_compositions, 3, 4)
