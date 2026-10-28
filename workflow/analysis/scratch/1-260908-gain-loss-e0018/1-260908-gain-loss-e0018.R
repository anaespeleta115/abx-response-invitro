source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260908-gain-loss-e0018/out/"


# this script takes most of the code from 1-260408-gainLossAME004.R and plots the loss in abundance of each family in the
# recipient community against the gain in abundance of each family in the mixture community.

# filtering steps are labeled throughout the script. These filters make the plots include or exclude taxa that colonized, but
# weren't originally present in the recipient community.

# ----------------- Compute Recipient Family Loss ------------------------

recipient_OTUs <- rbind(preabx_recipient_OTUs, postabx_recipient_OTUs)

recipient_fam_loss <- recipient_OTUs %>%
  group_by(mixture, Family, household) %>%
  summarise(
    fam_relAbundance = sum(relAbundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = mixture,
    values_from = fam_relAbundance,
    names_prefix = "recipient_fam_relAbundance_",
    values_fill = 0
  ) %>%
  # only keep families that were present pre-antibiotic
  filter(`recipient_fam_relAbundance_A-pre` > 0) %>%  # KEY FILTERING STEP # 1!!
  mutate(
    `recipient_fam_relAbundance_A-post_pseudo` = if_else(
    `recipient_fam_relAbundance_A-post` == 0,
    1e-4,
    `recipient_fam_relAbundance_A-post`
    ),
    recipient_ratio = `recipient_fam_relAbundance_A-pre` /
    `recipient_fam_relAbundance_A-post_pseudo`,
    recipient_log10_ratio = log10(recipient_ratio)
  )


# ----------------- Compute Post-Abx Mixture Family Gain ------------------------

fam_abundance_mixture <- mixture_OTUs %>%
  group_by(Family, household) %>%
  summarise(
    mixture_fam_abundance = sum(relAbundance),
    .groups = "drop"
  )


# get post-abx recipient abundance
Apost_recipient_fam_abundance <- recipient_OTUs %>%
  filter(
    mix_type == "A-post-1"
  ) %>%
  group_by(Family, household) %>%
  summarise(
    Apost_fam_abundance = sum(relAbundance),
    .groups = "drop"
  )


mixture_gain <- fam_abundance_mixture %>%
  left_join(
    Apost_recipient_fam_abundance,
    by = c("Family", "household")
  ) %>%
  mutate(
    Apost_fam_abundance = replace_na(Apost_fam_abundance, 1e-4),
    mix_ratio = mixture_fam_abundance / Apost_fam_abundance,
    mix_log10_ratio = log10(mix_ratio)
  )


# ----------------- Combine Gain and Loss Data ------------------------

compare_gain_loss <- mixture_gain %>%
  left_join(
    recipient_fam_loss,
    by = c("Family", "household")
  ) |> 
  mutate(
    recipient_ratio = replace_na(recipient_ratio, 1e-4),
    recipient_log10_ratio = log(recipient_ratio)
  ) |> # only keep families that were present pre-antibiotic
  filter(recipient_ratio != 1e-4)  # KEY FILTERING STEP # 2!!


# ----------------- Average Across Recipient's Mixes ------------------------


gain_loss_plot_e0018 <- compare_gain_loss %>%
  # filter(household %in% c("L", "P", "A", "Q")) |>
  ggplot(aes(x = recipient_log10_ratio, y = mix_log10_ratio, color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Recipient loss (log(A-pre recipient fam relAbundance / 
       A-post recipient fam relAbundance)", y = "Mixture gain (log(Mixture fam relAbundance / 
       A-post recipient fam relAbundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~household)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

gain_loss_plot_e0018

savePNGPDF(paste0(OUTDIR, "gain_loss_plot_e0018"), gain_loss_plot_e0018, 6, 6)




