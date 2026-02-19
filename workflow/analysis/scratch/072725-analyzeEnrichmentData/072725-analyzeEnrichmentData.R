# Questions I should be asking: 1.) Which families are showing up as enriched? Are they the same few families? Do the enriched OTUs 
# belong to a single family or are they all across the board? 

# 2.) Which subjects show the most enrichment? Why might that be? Can you relate that to the magnitude of species loss?

# 3.) How do those mixtures compare? What is the donor-recipient divergence?

# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")


# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072725-analyzeEnrichmentData/out/"

enrichment_summary_36 <- read.csv("C:/abx-response-invitro/data/enrichment_summary_day36.csv")

enrichment_ids <- read.csv("C:/abx-response-invitro/data/enrichment_summary_day36.csv")

enrichment_summary_64 <- read.csv("C:/abx-response-invitro/data/enrichment_summary_day64.csv")

# Set subject palette
SUBJECT <- c("XBA", "XDA", "XEA", "XKA")
PALETTE.SUBJECT <- c(  "#345995","#AD0505", "#daa520","#619E00")
names(PALETTE.SUBJECT) <- SUBJECT



# ---------------------------------------------- DAY 36 ----------------------------------------------------

# Prepare data
otu_enrichment_data <- enrichment_summary_36 %>%
  mutate(p_otu = -(p_value_otu)) %>%
  mutate(
    donor = str_split_fixed(mixture, "\\+", 2)[, 1],
    recipient = str_split_fixed(mixture, "\\+", 2)[, 2]
  ) %>%
  mutate(recipient = str_sub(recipient, 1, -5)) %>%
  mutate(donor = str_sub(donor, 1, -5)) %>% 
  select(mixture, donor, recipient, p_otu)

otu_enrichment_data <- otu_enrichment_data %>%
  group_by(donor) %>%
  mutate(mean_p_otu = mean(p_otu, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(donor = fct_reorder(donor, mean_p_otu))

otu_enrichment_plot <- otu_enrichment_data %>% ggplot(aes(x = donor, y = p_otu, color = recipient, group = recipient))+
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = PALETTE.SUBJECT)+
    labs(
      title = "",
      x = "Donor",
      y = expression(-log10(p[otu])),
      color = "Recipient"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      plot.title = element_text(size = 7, face = "bold")
    )+
    DEFAULTS.THEME_PRINT




savePNGPDF(paste0(OUTDIR, "enrichment_otu"), otu_enrichment_plot, 2, 3)



# Do the same for family enrichment

# FILTER OUT PRE-ABX MIXTURES

enrichment_summary_36 <- enrichment_summary_36 %>%
  mutate(day = str_sub(mixture, -3))
  # filter(day == "036")

# Prepare data
fam_enrichment_data <- enrichment_summary_36 %>%
  mutate(p_fam = -log10(p_value_fam)) %>%
  mutate(
    donor = str_split_fixed(mixture, "\\+", 2)[, 1],
    recipient = str_split_fixed(mixture, "\\+", 2)[, 2]
  ) %>%
  mutate(recipient = str_sub(recipient, 1, -5)) %>%
  mutate(donor = str_sub(donor, 1, -5)) %>% 
  select(mixture, donor, recipient, p_fam)

fam_enrichment_data <- fam_enrichment_data %>%
  group_by(donor) %>%
  mutate(mean_p_fam = mean(p_fam, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(donor = fct_reorder(donor, mean_p_fam))

fam_enrichment_plot <- fam_enrichment_data %>% ggplot(aes(x = donor, y = p_fam, color = recipient, group = recipient))+
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = PALETTE.SUBJECT)+
  labs(
    title = "",
    x = "Donor",
    y = expression(-log10(p[fam])),
    color = "Recipient"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  )+
  DEFAULTS.THEME_PRINT




savePNGPDF(paste0(OUTDIR, "enrichment_fam_36"), fam_enrichment_plot, 2, 3)


# ---------------------------------------------- DAY 64 ----------------------------------------------------



# Prepare data
otu_enrichment_data_64 <- enrichment_summary_64 %>%
  mutate(p_otu = -(p_value_otu)) %>%
  mutate(
    donor = str_split_fixed(mixture, "\\+", 2)[, 1],
    recipient = str_split_fixed(mixture, "\\+", 2)[, 2]
  ) %>%
  mutate(recipient = str_sub(recipient, 1, -5)) %>%
  mutate(donor = str_sub(donor, 1, -5)) %>% 
  select(mixture, donor, recipient, p_otu)

otu_enrichment_data_64 <- otu_enrichment_data_64 %>%
  group_by(donor) %>%
  mutate(mean_p_otu = mean(p_otu, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(donor = fct_reorder(donor, mean_p_otu))

otu_enrichment_plot_64 <- otu_enrichment_data_64 %>% ggplot(aes(x = donor, y = p_otu, color = recipient, group = recipient))+
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = PALETTE.SUBJECT)+
  labs(
    title = "",
    x = "Donor",
    y = expression(-log10(p[otu])),
    color = "Recipient"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  )+
  DEFAULTS.THEME_PRINT




savePNGPDF(paste0(OUTDIR, "enrichment_otu_64"), otu_enrichment_plot_64, 2, 3)



# Do the same for family enrichment

# Prepare data
fam_enrichment_data_64 <- enrichment_summary_64 %>%
  mutate(p_fam = -log10(p_value_fam)) %>%
  mutate(
    donor = str_split_fixed(mixture, "\\+", 2)[, 1],
    recipient = str_split_fixed(mixture, "\\+", 2)[, 2]
  ) %>%
  mutate(recipient = str_sub(recipient, 1, -5)) %>%
  mutate(donor = str_sub(donor, 1, -5)) %>% 
  select(mixture, donor, recipient, p_fam)

fam_enrichment_data_64 <- fam_enrichment_data_64 %>%
  group_by(donor) %>%
  mutate(mean_p_fam = mean(p_fam, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(donor = fct_reorder(donor, mean_p_fam))

fam_enrichment_plot_64 <- fam_enrichment_data_64 %>% ggplot(aes(x = donor, y = p_fam, color = recipient, group = recipient))+
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = PALETTE.SUBJECT)+
  labs(
    title = "",
    x = "Donor",
    y = expression(-log10(p[fam])),
    color = "Recipient"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  )+
  DEFAULTS.THEME_PRINT




savePNGPDF(paste0(OUTDIR, "enrichment_fam_64"), fam_enrichment_plot_64, 2, 3)


