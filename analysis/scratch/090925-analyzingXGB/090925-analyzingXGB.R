# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/090925-analyzingXDB/out/"


num_replicate <- 1

DAY <- c("029", "036")
PALETTE.DAY <- c("#FF68A1", "#8494FF")
names(PALETTE.DAY) <- DAY



# -------------------------- calculate species richness for XJB mixtures -----------------------------------------------


XJB_species_richness <- actual_colonizers_results %>% 
  filter(subject2 == "XJB", replicate == num_replicate, day == "036") %>% 
  group_by(mixture) %>% 
  summarize(species_richness = n_distinct(OTU))


# get colonization proportion for each mixture

XJB_colonization_prop <- colonization_prop_results %>% 
  filter(donor == "XJB", day != "064")


#------------------------ plot mixture colonization efficacy ------------------------------------------------------------



XJB_colonization_prop_plot <- XJB_colonization_prop %>% 
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


savePNGPDF(paste0(OUTDIR, "XJB_colonization_prop_plot_both"), XJB_colonization_prop_plot, 2, 3)



#------------------------ plot number of colonizers ------------------------------------------------------------

XJB_num_colonizers_plot <- XJB_colonization_prop %>% 
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


savePNGPDF(paste0(OUTDIR, "XJB_num_colonizers_plot_both"), XJB_num_colonizers_plot, 2, 3)

#------------------------ plot potential colonizers ------------------------------------------------------------



XJB_potential_colonizers_plot <- XJB_colonization_prop %>% 
  ggplot(aes(x = recipient, y = n_potential_colonizers, group = day, color = day)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    title = "",
    x = "Recipient",
    y = "Potential colonizers",
    color = "Day"
  ) +
  scale_y_continuous(limits = c(0, 110)) + 
  scale_color_manual(values = PALETTE.DAY)+
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "XJB_potential_colonizers_plot_both"), XJB_potential_colonizers_plot, 2, 3)


