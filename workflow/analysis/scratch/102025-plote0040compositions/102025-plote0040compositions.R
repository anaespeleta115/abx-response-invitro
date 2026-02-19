# load scripts
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102025-loade0040data/102025-loade0040data.R")

source("~/Documents/GitHub/abx-response-invitro/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102025-plote0040compositions/out/"



# Plot the total relative abundances (which are calculated per community) totaled up for each
p_compositions_1 <- ggplot(e0040 %>% filter(replicate == 1), aes(x = fct_relevel(community, "XEA-pre-abx", "XEA-post-abx-V1", "XEA-post-abx-V2"), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_compositions_2 <- ggplot(e0040 %>% filter(replicate == 2), aes(x = fct_relevel(community, "XEA-pre-abx", "XEA-post-abx-V1", "XEA-post-abx-V2"), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_compositions_3 <- ggplot(e0040 %>% filter(replicate == 3), aes(x = fct_relevel(community, "XEA-pre-abx", "XEA-post-abx-V1", "XEA-post-abx-V2"), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

p_compositions_4 <- ggplot(e0040 %>% filter(replicate == 4), aes(x = fct_relevel(community, "XEA-pre-abx", "XEA-post-abx-V1", "XEA-post-abx-V2"), y = relAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Community type", y = "", fill = "Family") +
  theme_minimal()+
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

combined_replicates <- p_compositions_1 / p_compositions_2 / p_compositions_3 / p_compositions_4

savePNGPDF(paste0(OUTDIR, "combined_replicates"), combined_replicates, 6, 3)

savePNGPDF(paste0(OUTDIR, "rep1"), p_compositions_1, 3, 3)

