source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-kxue-loadingScript/1-260406-kxue-loadingScript.R")

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

savePNGPDF(paste0(OUTDIR, "eACE010_community_compositions_byDay"), final_plot, 7, 7)

# Plot compositions faceted by their subject
eACE010_community_compositions <- eACE010_data %>% filter(round2plate == "A") %>% 
  ggplot(aes(x = factor(day), y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~subject)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eACE010_community_compositions

savePNGPDF(paste0(OUTDIR, "eACE010_community_compositions"), eACE010_community_compositions, 5, 5)



# Add richness information to the plots
richness <- data_eACE010_diversity %>% 
  rename(biosample = metadata) %>% 
  mutate(subject = str_sub(biosample, 1, 3),
         day = str_sub(biosample, -3),
         y_pos = 1.01) %>% 
  mutate(day = case_when(
    str_detect(day, "001") | str_detect(day, "002") | str_detect(day, "003") | str_detect(day, "022")| str_detect(day, "008") ~ "001",
    str_detect(day, "029") | str_detect(day, "028") | str_detect(day, "027") ~ "029",
    str_detect(day, "036") | str_detect(day, "037") | str_detect(day, "038") ~ "036",
    str_detect(day, "064")| str_detect(day, "063") | str_detect(day, "072") | str_detect(day, "059")| str_detect(day, "065") | str_detect(day, "066") ~ "064",
    TRUE ~ "0"
  ))
  

# eACE010_data_richness <- eACE010_data %>% 
#   left_join(data_eACE010_diversity, by = c("biosample"))

# Plot compositions faceted by their subject
eACE010_richness_compositions <- eACE010_data_richness %>% filter(round2plate == "A") %>% 
  ggplot(aes(x = factor(day), y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  geom_text(
    data = richness,
    aes(x = day, y = y_pos, label = nASVs),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_fill_manual(values = my_colors) +
  ylim(0, 1.1)+
  labs(x = "sample", y = "relative abundance")+
  facet_wrap(~subject)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

eACE010_richness_compositions

savePNGPDF(paste0(OUTDIR, "eACE010_richness_compositions"), eACE010_richness_compositions, 5, 5)


