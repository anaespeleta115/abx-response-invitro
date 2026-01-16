source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/011426gete0029LostTaxa/out/"


XBA_ASVs <- recipient_ASVs %>% 
  filter(subject1 == "XKA")

e0029_day29 <- unique(XBA_ASVs %>% 
                          filter(replicate == 1, day == "029") %>% 
                          pull(OTU))

e0029_day36 <- unique(XBA_ASVs %>% 
                             filter(replicate == 1, day == "036") %>% 
                             pull(OTU))


# post_abx_lost_taxa_v1 <- setdiff(e0040_pre_abx, e0040_postabx_v1)
# post_abx_lost_taxa_v2 <- setdiff(e0040_pre_abx, e0040_postabx_v2)


post_abx_lost_taxa <- setdiff(e0029_day29, e0029_day36) # get taxa present in pre-abx community and not in post-abx community

e0029_recipients_lost <-recipient_ASVs %>% filter(OTU %in% post_abx_lost_taxa) %>%
  select(OTU, Family) %>%
  distinct()


# --------------------------------- PLOT DAY 36 LOST TAXA ---------------------------------------
# For coloring only the lost taxa
## alpha = factor(lost)

p_pre_abx_lost <- ggplot(recipient_ASVs %>% filter(replicate == 1, day == "029", subject1 == "XKA") %>% mutate(lost = ifelse(OTU %in% post_abx_lost_taxa, 1, 0)), aes(x = subject1, y = relAbundance, fill = Family,  alpha = factor(lost))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  labs(title = "" , x = "", y = "", fill = "Family") +
  theme(
    legend.position = "none"
  ) +
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "compositions_pre-abx_XKA_lost"), p_pre_abx_lost, 2, 1)


