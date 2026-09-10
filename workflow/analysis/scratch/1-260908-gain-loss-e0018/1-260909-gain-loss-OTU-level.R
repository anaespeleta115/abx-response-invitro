source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260908-gain-loss-e0018/out/"


# this script adjusts the gain/loss analysis from 1-260908-gain-loss-e0018 and changes it to compute gain and loss not at
# the family level, but at the OTU level. This will hopefully tell us whether the patterns we are observing are only occuring at
# the family level and not at the OTU level. 



# ----------------- Compute Recipient Family Loss ------------------------

recipient_OTUs <- rbind(preabx_recipient_OTUs, postabx_recipient_OTUs)

recipient_OTU_loss <- recipient_OTUs %>%
  group_by(Family, OTU, household, mix_type) %>%
  summarise(
    OTU_relAbundance = sum(relAbundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = mix_type,
    values_from = OTU_relAbundance,
    names_prefix = "recipient_OTU_relAbundance_",
    values_fill = 0
  ) %>%
  # only keep families that were present pre-antibiotic
  filter(`recipient_OTU_relAbundance_A-pre-1` > 0) %>%  # KEY FILTERING STEP !!
  mutate(
    `recipient_OTU_relAbundance_A-post-1_pseudo` = if_else(
      `recipient_OTU_relAbundance_A-post-1` == 0,
      1e-4,
      `recipient_OTU_relAbundance_A-post-1`
    ),
    recipient_ratio = `recipient_OTU_relAbundance_A-pre-1` /
      `recipient_OTU_relAbundance_A-post-1_pseudo`,
    recipient_log10_ratio = log10(recipient_ratio)
  )


# ----------------- Compute Post-Abx Mixture Family Gain ------------------------

OTU_abundance_mixture <- mixture_OTUs %>%
  group_by(Family,OTU, household) %>%
  summarise(
    mixture_OTU_abundance = sum(relAbundance),
    .groups = "drop"
  )


# get post-abx recipient abundance
Apost_recipient_OTU_abundance <- recipient_OTUs %>%
  filter(
    mix_type == "A-post-1"
  ) %>%
  group_by(Family, OTU, household) %>%
  summarise(
    Apost_OTU_abundance = sum(relAbundance),
    .groups = "drop"
  )


mixture_gain <- OTU_abundance_mixture %>%
  left_join(
    Apost_recipient_OTU_abundance,
    by = c("OTU", "household")
  ) %>%
  mutate(
    Apost_OTU_abundance = replace_na(Apost_OTU_abundance, 1e-4),
    mix_ratio = mixture_OTU_abundance / Apost_OTU_abundance,
    mix_log10_ratio = log10(mix_ratio)
  )


# ----------------- Combine Gain and Loss Data ------------------------

compare_gain_loss <- mixture_gain %>%
  left_join(
    recipient_OTU_loss,
    by = c("OTU", "household")
  ) 

compare_gain_loss <- compare_gain_loss |>  
  mutate(
    `recipient_OTU_relAbundance_A-pre-1` = replace_na(`recipient_OTU_relAbundance_A-pre-1`, 1e-4),
    `recipient_OTU_relAbundance_A-post-1` = replace_na(`recipient_OTU_relAbundance_A-post-1`, 1e-4),
    `recipient_OTU_relAbundance_A-post-1_pseudo` = replace_na(`recipient_OTU_relAbundance_A-post-1_pseudo`, 1e-4),
    recipient_ratio = replace_na(recipient_ratio, 1e-4),
    recipient_log10_ratio = log(recipient_ratio)
  ) # make sure to convert all NAs (taxa not present in recipient) with pseudocount


# ----------------- Average Across Recipient's Mixes ------------------------


otu_gain_loss_plot_e0018 <- compare_gain_loss %>%
  ggplot(aes(x = recipient_log10_ratio, y = mix_log10_ratio, color = Family.x))+
  geom_point(size = 0.5, alpha = 1)+
  # geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Recipient loss (log(A-pre recipient OTU relAbundance / 
       A-post recipient OTU relAbundance)", y = "Mixture gain (log(Mixture OTU relAbundance / 
       A-post recipient OTU relAbundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~household)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

otu_gain_loss_plot_e0018

savePNGPDF(paste0(OUTDIR, "otu_gain_loss_plot_e0018"), otu_gain_loss_plot_e0018, 6, 6)


# ----------------- Correlate Recipient 29 and Mixture 36 at the OTU-level -----------



# ----------------- Compute Pre-abx Recipient Abundance (recipient 29) ------------------------

recipient_OTU_abundance_preabx <- preabx_recipient_OTUs %>%
  group_by(Family, OTU, household) %>%
  summarise(recipient_OTU_abundance = sum(relAbundance))

# ----------------- Compute Post-Abx Family Abundance (mixture 36) -------------------------------

OTU_abundance_mixture <- mixture_OTUs %>%
  group_by(OTU, household) %>%
  summarise(
    mixture_OTU_abundance = sum(relAbundance),
    .groups = "drop"
  )

# ----------------- Put both family abundance datasets together ------------------------

compare_recipient_mixture <- recipient_OTU_abundance_preabx %>% 
  left_join(OTU_abundance_mixture, by = c("household", "OTU"), relationship = "many-to-many") %>% 
  replace_na(list(recipient_OTU_abundance = 1e-4, mixture_OTU_abundance = 1e-4))


# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

recipient_mixture_plot <- compare_recipient_mixture %>% 
  # filter(household %in% c("L", "P", "A", "Q")) |> 
  ggplot(aes(x = recipient_OTU_abundance, y = mixture_OTU_abundance, color = Family))+
  geom_point(size = 0.5, alpha = 1)+
  # geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "log(Pre-abx recipient OTU abundance)", y = "log(Post-abx mixture OTU abundance)")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~household)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

recipient_mixture_plot

savePNGPDF(paste0(OUTDIR, "otu_recipient_mixture_plot"), recipient_mixture_plot, 6, 6)

