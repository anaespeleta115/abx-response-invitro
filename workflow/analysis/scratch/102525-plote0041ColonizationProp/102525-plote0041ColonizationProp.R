# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-gete0041colonization/102125-gete0041colonization.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102525-plote0041colonizationProp/out/"



VERSION <- c("pre-abx", "post-abx-V2", "post-abx-V1")
PALETTE.VERSION <- c(paletteer_d("nationalparkcolors::ArcticGates"))
names(PALETTE.VERSION) <- VERSION


mixture_colonization_data <- mixture_colonization_full %>% 
  filter(donor != "XBB-029" & donor != "XEB-029" & donor != "XDB-029") %>% 
  group_by(mixture) %>% 
  summarize(total_colonizers = sum(colonizer)) %>% 
  left_join(mixture_colonization_full %>% select(mixture, recipient, donor), by = "mixture") %>% 
  distinct()


colonization_plot <-  ggplot(mixture_colonization_data, aes(x = fct_relevel(recipient, c("pre-abx", "post-abx-V1", "post-abx-V2")), y = total_colonizers, fill = recipient)) +
  geom_boxplot() +
  labs(
    title = "",
    x = "Recipient",
    y = "Number of colonizers",
    fill = "Pre/Post-abx"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = PALETTE.SUBJECT)+
  DEFAULTS.THEME_PRINT+
  theme(legend.position = "none")

savePNGPDF(paste0(OUTDIR, "mixtureColonization_adjusted"), colonization_plot, 3, 3)


prop_colonization_plot <- 
  ggplot(mixture_colonization_proportion_full %>% filter(donor != "XBB-029" & donor != "XEB-029" & donor != "XDB-029"), aes(x = fct_relevel(recipient, c("pre-abx", "post-abx-V1", "post-abx-V2")), y = colonization_proportion, fill = recipient)) +
  geom_boxplot() +
  labs(
    title = "",
    x = "Recipient",
    y = "Colonization proportion",
    fill = "Pre/Post-abx"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = PALETTE.VERSION)+
  DEFAULTS.THEME_PRINT+
  theme(legend.position = "none")

savePNGPDF(paste0(OUTDIR, "mixtureColonizationProp_adjusted"), prop_colonization_plot, 2, 2)