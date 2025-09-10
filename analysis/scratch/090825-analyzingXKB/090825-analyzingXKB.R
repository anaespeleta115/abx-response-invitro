# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/090825-analyzingXKB/out/"


num_replicate <- 1

DAY <- c("029", "036")
PALETTE.DAY <- c("#00BA38", "#F37735")
names(PALETTE.DAY) <- DAY



# -------------------------- calculate species richness for XKB mixtures -----------------------------------------------


XKB_species_richness <- actual_colonizers_results %>% 
  filter(subject2 == "XKB", replicate == num_replicate, day == "036") %>% 
  group_by(mixture) %>% 
  summarize(species_richness = n_distinct(OTU))


# get colonization proportion for each mixture

XKB_colonization_prop <- colonization_prop_results %>% 
  filter(donor == "XKB", day != "064")


#------------------------ plot mixture colonization efficacy ------------------------------------------------------------



XKB_colonization_prop_plot <- XKB_colonization_prop %>% 
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
  scale_color_manual(values = PALETTE.DAY)+
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XKB_colonization_prop_plot_both"), XKB_colonization_prop_plot, 2, 3)



#------------------------ plot number of colonizers ------------------------------------------------------------

XKB_num_colonizers_plot <- XKB_colonization_prop %>% 
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
  scale_color_manual(values = PALETTE.DAY)+ 
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XKB_num_colonizers_plot_both"), XKB_num_colonizers_plot, 2, 3)

#------------------------ plot potential colonizers ------------------------------------------------------------



XKB_potential_colonizers_plot <- XKB_colonization_prop %>% 
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
  scale_color_manual(values = PALETTE.DAY)+
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XKB_potential_colonizers_plot_both"), XKB_potential_colonizers_plot, 2, 3)


