# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")




subject_composition <- function(mix_id, donor_id, recipient_id, replicate) {
  
  # Subset mixture and recipient ASVs for the given replicate
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == replicate)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == replicate)
  
  mixture_ASVs_mix <- actual_colonizers_results %>% 
    filter(mixture == mix_id, replicate == 1) %>% 
    ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Mixture
      ",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  
  mixture_ASVs_colonized <- actual_colonizers_results %>% 
    filter(mixture == mix_id, replicate == 1) %>% 
    ggplot(aes(x = mixture, y = relAbundance, fill = Family, alpha = factor(actual_colonizer))) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1)) +  # colonizer == 0 → transparent
    labs(
      title = "Actual 
Colonizers",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  mixture_ASVs_differential <- actual_colonizers_results %>%
    filter(mixture == mix_id, replicate == replicate, day == "064") %>%
    ggplot(aes(x = mixture, y = relAbundance, fill = Family, alpha = factor(diff_colonizer_36))) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1)) +  # colonizer == 0 → transparent
    labs(
      title = "Differential
Colonizers",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(),
          axis.text.x = element_text(size = 10))
  
  
  lost_ASVs <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == 1, day == "029") %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family, alpha = factor(lost_strain_29_64)))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
    labs(
      title = "Lost 
Strains",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  recipient_ASVs <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == 1) %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Recipient
      ",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  donor_ASVs <- single_donor_ASVs %>% 
    filter(biosample1 == donor_id)%>%
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Donor
      ",
      x = "",
      y = "Total Relative Abundance"
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  
  row_of_plots <- plot_grid(
    donor_ASVs,
    recipient_ASVs,
    mixture_ASVs_mix,
    mixture_ASVs_colonized,
    lost_ASVs,
    mixture_ASVs_differential,
    ncol = 6
  )
  
  # Add a centered title above the row
  final_plot <- plot_grid(
    ggdraw() + draw_label("", fontface = "bold", size = 12, hjust = 0.5),
    row_of_plots,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  return(final_plot)
}


# Function call
subject_composition_plot_day29 <- subject_composition("XDB-029+XEA-029", "XDB-029", "XEA-029", 1) + DEFAULTS.THEME_PRINT
subject_composition_plot_day64 <- subject_composition("XDB-029+XEA-064", "XDB-029", "XEA-064", 1) + DEFAULTS.THEME_PRINT


combined_plot <- subject_composition_plot_day29 /
  subject_composition_plot_day64



