# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-getColonization/072125-getColonization.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")

subject_composition <- function(mix_id, donor_id, recipient_id, replicate) {
  
  # Subset mixture and recipient ASVs for the given replicate
  mix_asvs_subset <- mixture_ASVs %>% 
    filter(mixture == mix_id, replicate == replicate)
  
  recipient_asvs_subset <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == replicate)
  
  mixture_ASVs_mix <- colonization_results %>% 
    filter(mixture == mix_id, replicate == replicate) %>% 
    ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Mixture",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 7))
  
  
  mixture_ASVs_colonized <- colonization_results %>% 
    filter(mixture == mix_id, replicate == replicate) %>% 
    ggplot(aes(x = mixture, y = relAbundance, fill = Family, alpha = factor(colonization))) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1)) +  # colonizer == 0 → transparent
    labs(
      title = "Colonized",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 7))
  
  
  recipient_ASVs <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id, replicate == replicate)%>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Recipient",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 7))
  
  donor_ASVs <- single_donor_ASVs %>% 
    filter(biosample1 == donor_id)%>%
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Donor",
      x = "",
      y = "Total Relative Abundance"
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 7))
  
  
  row_of_plots <- plot_grid(
    donor_ASVs,
    recipient_ASVs,
    mixture_ASVs_mix,
    mixture_ASVs_colonized,
    ncol = 4
  )
  
  # Add a centered title above the row
  final_plot <- plot_grid(
    ggdraw() + draw_label("", fontface = "bold", size = 12, hjust = 0.5),
    row_of_plots,
    ncol = 1,
    rel_heights = c(0.1, 1)  # Reserve top space for title
  )
  
  return(final_plot)
}


# Function calls
subject_composition_plot_day29 <- subject_composition("XKB-029+XEA-029", "XKB-029", "XEA-029", 1) + DEFAULTS.THEME_PRINT
subject_composition_plot_day36 <- subject_composition("XKB-029+XEA-036", "XKB-029", "XEA-036", 1) + DEFAULTS.THEME_PRINT
subject_composition_plot_day64 <- subject_composition("XKB-029+XEA-064", "XKB-029", "XEA-064", 1) + DEFAULTS.THEME_PRINT

# Combined day plots
combined_plot <- subject_composition_plot_day29 /
  subject_composition_plot_day36 /
  subject_composition_plot_day64