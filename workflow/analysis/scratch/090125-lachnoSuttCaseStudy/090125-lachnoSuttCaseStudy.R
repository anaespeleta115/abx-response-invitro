# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/090125-lachnoSuttCaseStudy/out/"


num_replicate <- 1



recipient_composition <- function(recipient_id) {
  
  
  recipient_ASVs_29 <- recipient_ASVs %>% 
    mutate(day = str_sub(biosample1, -3),
           subject1 = str_sub(biosample1, 1, -5)) %>% 
    filter(day == "029", subject1 == recipient_id, replicate == num_replicate) %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family, alpha = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(
      values = c("Lachnospiraceae" = 1, 
                 "Sutterellaceae" = 1, 
                 "Bacteroidaceae" = 1),
      na.value = 0
    ) +
    labs(
      title = "",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  recipient_ASVs_36 <- recipient_ASVs %>%
    mutate(day = str_sub(biosample1, -3),
           subject1 = str_sub(biosample1, 1, -5)) %>% 
    filter(day == "036", subject1 == recipient_id, replicate == num_replicate) %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family, alpha = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
    scale_alpha_manual(
      values = c("Lachnospiraceae" = 1, 
                 "Sutterellaceae" = 1, 
                 "Bacteroidaceae" = 1),
      na.value = 0
    ) +
    labs(
      title = "",
      x = "",
      y = ""
    ) +
    theme(legend.position = "none",
          axis.line=element_blank(), 
          axis.text.x = element_text(size = 10))
  
  
  recipient_plot <- plot_grid(
    recipient_ASVs_29,
    recipient_ASVs_36,
    ncol = 2
  )
  
  # Make a list of mixture plots, one per biosample2
  mixture_plots <- actual_colonizers_results %>%
    filter(day == "036", subject1 == recipient_id, replicate == num_replicate) %>%
    split(.$biosample2) %>%
    map(~ ggplot(.x, aes(x = mixture, y = relAbundance, fill = Family, alpha = Family)) +
          geom_bar(stat = "identity", color = "black") +
          scale_fill_manual(values = my_colors) +
          scale_alpha_manual(
            values = c("Lachnospiraceae" = 1, 
                       "Sutterellaceae" = 1, 
                       "Bacteroidaceae" = 1),
            na.value = 0
          ) +
          labs(x = "", y = "", title = unique(.x$biosample2)) +
          theme(
            legend.position = "none",
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          ))
  
  
  # Stitch them into one row
  mixture_plot <- plot_grid(plotlist = mixture_plots, nrow = 1)
  
  # Add a centered title above the row
  final_plot <- plot_grid(
    recipient_plot,
    mixture_plot,
    nrow = 1,
    rel_widths = c(1, 5)  # make Mixture twice as wide
  )
  
  return(final_plot)
}


recipient_mixtures_36_plot <- recipient_composition("XKA")

savePNGPDF(paste0(OUTDIR, "XKA_lachnoBacSutt_plot"), recipient_mixtures_36_plot, 4, 20)

