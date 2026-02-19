# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/090125-analyzingXHB/out/"


num_replicate <- 1


# -------------------------- calculate species richness for XHB mixtures -----------------------------------------------


XHB_species_richness <- actual_colonizers_results %>% 
  filter(subject2 == "XHB", replicate == num_replicate, day == "036") %>% 
  group_by(mixture) %>% 
  summarize(species_richness = n_distinct(OTU))


# get colonization proportion for each mixture


XHB_colonization_prop <- colonization_prop_results %>% 
  filter(donor == "XHB", day != "064")
# 
# XHB_colonization_prop_36 <- colonization_prop_results %>% 
#   filter(donor == "XHB", day == "036")



#------------------------ plot mixture colonization efficacy ------------------------------------------------------------




XHB_colonization_prop_plot <- XHB_colonization_prop %>% 
  ggplot(aes(x = recipient, y = prop_colonizers, group = day, color = day)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
      title = "",
      x = "Recipient",
      y = "Colonization efficacy",
      color = "Day"
  ) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XHB_colonization_prop_plot_both"), XHB_colonization_prop_plot, 2, 3)



#------------------------ plot number of colonizers ------------------------------------------------------------

XHB_num_colonizers_plot <- XHB_colonization_prop %>% 
  ggplot(aes(x = recipient, y = n_actual_colonizers, group = day, color = day)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    title = "",
    x = "Recipient",
    y = "Num colonizers",
    color = "Day"
  ) +
  scale_y_continuous(limits = c(0, 100)) + 
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XHB_num_colonizers_plot_both"), XHB_num_colonizers_plot, 2, 3)

#------------------------ plot potential colonizers ------------------------------------------------------------



XHB_potential_colonizers_plot <- XHB_colonization_prop %>% 
  ggplot(aes(x = recipient, y = n_potential_colonizers, group = day, color = day)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    title = "",
    x = "Recipient",
    y = "Potential colonizers",
    color = "Day"
  ) +
  scale_y_continuous(limits = c(0, 100)) + 
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XHB_potential_colonizers_plot_both"), XHB_potential_colonizers_plot, 2, 3)