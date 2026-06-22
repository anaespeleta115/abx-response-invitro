# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004colonization/1-260428-geteAME004colonization.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260408-gainLossAME004/out/"

curr_replicate <- 1


# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_loss <- recipient_ASVs %>%
  filter(
    day %in% c("029", "036"),
    replicate == curr_replicate,
    relAbundance > 1e-3
  ) %>%
  group_by(day, Family, subject) %>%
  summarise(
    fam_relAbundance = sum(relAbundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = day,
    values_from = fam_relAbundance,
    names_prefix = "recipient_fam_relAbundance_day_",
    values_fill = 0
  ) %>%
  # keep families that were present pre-antibiotic
  filter(recipient_fam_relAbundance_day_029 > 0) %>%
  mutate(
    recipient_fam_relAbundance_day_036_pc = if_else(
      recipient_fam_relAbundance_day_036 == 0,
      1e-4,
      recipient_fam_relAbundance_day_036
    ),
    recipient_ratio = recipient_fam_relAbundance_day_029 /
      recipient_fam_relAbundance_day_036_pc,
    recipient_log10_ratio = log10(recipient_ratio)
  )


# ----------------- Compute Post-Abx Family Gain ------------------------

fam_abundance_mixture <- colonization_results %>%
  filter(
    day == "036",
    replicate == curr_replicate,
    relAbundance > 1e-3
  ) %>%
  mutate(subject2 = str_sub(biosample2, 1, 3)) %>%
  group_by(Family, subject, subject2) %>%
  summarise(
    mixture_abundance36 = sum(relAbundance),
    .groups = "drop"
  )


recipient_fam_36 <- recipient_ASVs %>%
  filter(
    day == "036",
    replicate == curr_replicate,
    relAbundance > 1e-3
  ) %>%
  group_by(Family, subject) %>%
  summarise(
    recipient_abundance36 = sum(relAbundance),
    .groups = "drop"
  )


gain_mixture <- fam_abundance_mixture %>%
  left_join(
    recipient_fam_36,
    by = c("Family", "subject")
  ) %>%
  mutate(
    recipient_abundance36 = replace_na(recipient_abundance36, 1e-4),
    mix_ratio = mixture_abundance36 / recipient_abundance36,
    mix_log10_ratio = log10(mix_ratio)
  )


# ----------------- Combine Gain and Loss Data ------------------------

compare_gain_loss <- gain_mixture %>%
  left_join(
    recipient_fam_loss,
    by = c("Family", "subject")
  )


# ----------------- Average Across Recipient's Mixes ------------------------

gain_loss_avg <- gain_mixture %>%
  group_by(Family, subject) %>%
  summarise(
    avg_mix_log10_ratio = mean(mix_log10_ratio, na.rm = TRUE),
    n_mixes = n(),
    .groups = "drop"
  ) %>%
  left_join(
    recipient_fam_loss,
    by = c("Family", "subject")
  ) 
  ## Uncomment this to display taxa that colonized, though not present in day 29 recipient
  # mutate(
  #   recipient_log10_ratio = replace_na(recipient_log10_ratio, 0),
  #   recipient_fam_relAbundance_day_029 = replace_na(recipient_fam_relAbundance_day_029, 0)
  # )



gain_loss_avg_plot <- gain_loss_avg %>%
  ggplot(aes(x = recipient_log10_ratio, y = avg_mix_log10_ratio, color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "Recipient loss (log(recipient day 29 relAbundance / 
       recipient day 36 relAbundance)", y = "Average mixture gain (log(mixture day 36 relAbundance / 
       recipient 36 relAbundance")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

gain_loss_avg_plot

savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot"), gain_loss_avg_plot, 3, 4)



# ---------------- Get the unexpected colonizers -----------------

unexpected_colonizers <- gain_loss_avg %>% 
  filter(is.na(recipient_fam_relAbundance_day_029)) %>% 
  group_by(subject, Family) %>% 
  summarize(total_relAbundance = sum(avg_mix_log10_ratio))

unexpected_colonizers_prev <- unexpected_colonizers %>% 
  group_by(Family, total_relAbundance) %>% 
  summarize(num_subjects = n(), avg_abundance = mean(total_relAbundance))


