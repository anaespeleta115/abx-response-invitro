# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/090925-analyzingXGB/out/"


num_replicate <- 1

DAY <- c("029", "036")
PALETTE.DAY <- c("#FF68A1", "#8494FF")
names(PALETTE.DAY) <- DAY



# -------------------------- calculate species richness for XCB mixtures -----------------------------------------------


XCB_species_richness <- actual_colonizers_results %>% 
  filter(subject2 == "XCB", replicate == num_replicate, day == "036") %>% 
  group_by(mixture) %>% 
  summarize(species_richness = n_distinct(OTU))


# get colonization proportion for each mixture

XCB_colonization_prop <- colonization_prop_results %>% 
  filter(donor == "XCB", day != "064")


#------------------------ plot mixture colonization efficacy ------------------------------------------------------------



XCB_colonization_prop_plot <- XCB_colonization_prop %>% 
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


savePNGPDF(paste0(OUTDIR, "XCB_colonization_prop_plot_both"), XCB_colonization_prop_plot, 2, 3)



#------------------------ plot number of colonizers ------------------------------------------------------------

XCB_num_colonizers_plot <- XCB_colonization_prop %>% 
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


savePNGPDF(paste0(OUTDIR, "XCB_num_colonizers_plot_both"), XCB_num_colonizers_plot, 2, 3)

#------------------------ plot potential colonizers ------------------------------------------------------------



XCB_potential_colonizers_plot <- XCB_colonization_prop %>% 
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


savePNGPDF(paste0(OUTDIR, "XCB_potential_colonizers_plot_both"), XCB_potential_colonizers_plot, 2, 3)


