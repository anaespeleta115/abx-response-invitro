# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/092525-compareLossGainDifference/out/"

curr_replicate <- 1

library(ppcor)

# ----------------- Compute Recipient Family Loss ------------------------

recipient_fam_loss <- recipient_ASVs %>%
  filter(day == "029" | day == "036", replicate == curr_replicate, relAbundance > 1e-4) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(day, Family, subject1) %>%
  summarise(fam_relAbundance = sum(relAbundance))

recipient_fam_loss <- recipient_fam_loss %>%
  pivot_wider(
    names_from = day,
    values_from = fam_relAbundance,
    names_prefix = "recipient_fam_relAbundance_day_"
  )%>%
  replace_na(list(recipient_fam_relAbundance_day_029 = 0, recipient_fam_relAbundance_day_036 = 0)) %>%
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(recipient_fam_relAbundance_day_029 == 0, 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(recipient_fam_relAbundance_day_036 == 0, 1e-4, recipient_fam_relAbundance_day_036),
    # recipient_ratio = recipient_fam_relAbundance_day_029 / recipient_fam_relAbundance_day_036   # check this !!!
    recipient_diff = log10(recipient_fam_relAbundance_day_029)
  )

#X = log(recipient_29)-log(recipient_36)

#Y = log(mixture_36)-log(mixture_36)

# ----------------- Compute Post-Abx Family Gain (recipient 36 -> mixture 36) -------------------------------

fam_abundance_mixture <- actual_colonizers_results %>% 
  filter(day == "036", replicate == curr_replicate) %>% 
  group_by(day, Family, subject1, subject2) %>% 
  summarise(fam_relAbundance = sum(relAbundance)) 

recipient_fam_36 <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>%
  group_by(day, Family, subject1) %>%
  summarise(fam_relAbundance = sum(relAbundance)) 

gain_mixture <- fam_abundance_mixture %>%
  left_join(recipient_fam_36, by = c("subject1", "Family"), relationship = "many-to-many") %>% 
  replace_na(list(fam_relAbundance.y = 0, fam_relAbundance.x = 0)) %>%
  mutate(
    fam_relAbundance.y = if_else(fam_relAbundance.y == 0, 1e-4, fam_relAbundance.y),
    fam_relAbundance.x = if_else(fam_relAbundance.x == 0, 1e-4, fam_relAbundance.x),
    # mix_ratio = fam_relAbundance.x / fam_relAbundance.y # x = mix 36, y = recipient 36
    mix_diff = log10(fam_relAbundance.x)
  ) 


# ----------------- Combine gain and loss data --------------------------------------------------------------

compare_gain_loss <- gain_mixture %>% 
  left_join(recipient_fam_loss, by = c("Family", "subject1")) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(is.na(recipient_fam_relAbundance_day_029), 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(is.na(recipient_fam_relAbundance_day_036), 1e-4, recipient_fam_relAbundance_day_036),
    recipient_diff = if_else(is.na(recipient_diff), 1e-4, recipient_diff)
  ) %>% 
  filter(recipient_fam_relAbundance_day_029 != 1e-04 & recipient_fam_relAbundance_day_036 != 1e-04)


# ----------------- Average across a recipient's mixes fold change for each recipient ------------------------

gain_loss_avg <- gain_mixture %>% 
  group_by(Family, subject1) %>% 
  summarise(avg_mix_diff = mean(mix_diff)) %>% 
  left_join(recipient_fam_loss, by = c("Family", "subject1")) %>% 
  mutate(
    recipient_fam_relAbundance_day_029 = if_else(is.na(recipient_fam_relAbundance_day_029), 1e-4, recipient_fam_relAbundance_day_029),
    recipient_fam_relAbundance_day_036 = if_else(is.na(recipient_fam_relAbundance_day_036), 1e-4, recipient_fam_relAbundance_day_036),
    recipient_diff = if_else(is.na(recipient_diff), 1e-4, recipient_diff)
  )
  # filter(recipient_diff != 1e-04)

gain_loss_avg_plot <- gain_loss_avg %>% 
  ggplot(aes(x = recipient_diff, y = avg_mix_diff, color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  # geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) + # correlation line
  labs(x = "log(recipient 29 fam relAbundance)", y = "log(mixture 36 fam relAbundance)")+
  scale_color_manual(values = my_colors)+
  facet_wrap(~subject1)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "gain_loss_avg_plot_numerators_log"), gain_loss_avg_plot, 3, 4)


# ----------------- Run stats ------------------------

compare_gain_loss <- as.data.frame(compare_gain_loss)
compare_gain_loss <- compare_gain_loss %>%
  mutate(
    eps = 1e-6,
    U = log(recipient_fam_relAbundance_day_029 + eps),
    V = log(fam_relAbundance.x + eps),
    W = log(recipient_fam_relAbundance_day_036 + eps),
    X = U - W,   # current x
    Y = V - W,   # current y
    D = V - U    # difference of numerators; no shared denominator
  )

coef(lm(Y ~ X, compare_gain_loss))         # slope ~ 1 suggests denominator-driven
cor(compare_gain_loss$U, compare_gain_loss$V, method="spear")  # relation of numerators only
# partial correlation of V and U controlling for W:
ppcor::pcor.test(compare_gain_loss$V, compare_gain_loss$U, compare_gain_loss$W, method="spear")


