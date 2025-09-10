# Load data
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/083125-plotRecipientCompositions/out/"


num_replicate <- 1



recipient_composition <- function(recipient_id) {

  
  recipient_ASVs_29 <- recipient_ASVs %>% 
    mutate(day = str_sub(biosample1, -3),
           subject1 = str_sub(biosample1, 1, -5)) %>% 
    filter(day == "029", subject1 == recipient_id, replicate == num_replicate) %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
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
    ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = my_colors) +
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
    map(~ ggplot(.x, aes(x = mixture, y = relAbundance, fill = Family)) +
          geom_bar(stat = "identity", color = "black") +
          scale_fill_manual(values = my_colors) +
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

savePNGPDF(paste0(OUTDIR, "XKA_recipient_mixtures_36_plot"), recipient_mixtures_36_plot, 4, 20)


# ----------------------------------------PLOT XHB-----------------------------------------------------


XHB_mixture_plots <- actual_colonizers_results %>%
  filter(day == "036", subject2 == "XHB", replicate == num_replicate) %>%
  split(.$biosample1) %>%
  map(~ggplot(.x, aes(x = mixture, y = relAbundance, fill = Family, alpha = factor(diff_colonizer_36))) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = my_colors)  +
        scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
        labs(x = "", y = "", title = "") +
        facet_wrap(~biosample1)+
        theme(
          legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ))

# Stitch them into one row
XHB_mixture_plots <- plot_grid(plotlist = XHB_mixture_plots, nrow = 1)

XHB <- single_donor_ASVs %>%
  mutate(day = str_sub(biosample1, -3),
         subject1 = str_sub(biosample1, 1, -5)) %>% 
  filter(subject1 == "XHB", replicate == num_replicate) %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.text.x = element_text(size = 10))


XHB_plots <- plot_grid(
  XHB,
  XHB_mixture_plots,
  nrow = 1,
  rel_widths = c(1, 3)  # make Mixture twice as wide
)


savePNGPDF(paste0(OUTDIR, "XHB_diff_colonization_plots_36"), XHB_plots, 4, 10)



# ----------------------------------------PLOT XKB-----------------------------------------------------


XKB_mixture_plots <- actual_colonizers_results %>%
  filter(day == "036", subject2 == "XKB", replicate == num_replicate) %>%
  split(.$biosample1) %>%
  map(~ggplot(.x, aes(x = mixture, y = relAbundance, fill = Family, alpha = factor(actual_colonizer))) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = my_colors) +
        scale_alpha_manual(values = c("0" = 0, "1" = 1))+
        labs(x = "", y = "", title = "") +
        facet_wrap(~biosample1)+
        theme(
          legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ))

# Stitch them into one row
XKB_mixture_plots <- plot_grid(plotlist = XKB_mixture_plots, nrow = 1)

XKB <- single_donor_ASVs %>%
  mutate(day = str_sub(biosample1, -3),
         subject1 = str_sub(biosample1, 1, -5)) %>% 
  filter(subject1 == "XKB", replicate == num_replicate) %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.text.x = element_text(size = 10))


XKB_plots <- plot_grid(
  XKB,
  XKB_mixture_plots,
  nrow = 1,
  rel_widths = c(1, 3)  # make Mixture twice as wide
)


savePNGPDF(paste0(OUTDIR, "XKB_plots_colonization_36"), XKB_plots, 4, 10)


# ----------------------------------------PLOT ALL DONORS -----------------------------------------------------

donor <- "XJB"

# ADD FOR COLONIZATION
## alpha = factor(actual_colonizer)

XJB_mixture_plots <- actual_colonizers_results %>%
  filter(day == "036", subject2 == donor, replicate == num_replicate) %>%
  split(.$biosample1) %>%
  map(~ggplot(.x, aes(x = mixture, y = relAbundance, fill = Family)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = my_colors) +
        # scale_alpha_manual(values = c("0" = 0, "1" = 1))+
        labs(x = "", y = "", title = "") +
        facet_wrap(~biosample1)+
        theme(
          legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ))

# Stitch them into one row
XJB_mixture_plots <- plot_grid(plotlist = XJB_mixture_plots, nrow = 1)

XJB <- single_donor_ASVs %>%
  mutate(day = str_sub(biosample1, -3),
         subject1 = str_sub(biosample1, 1, -5)) %>% 
  filter(subject1 == donor, replicate == num_replicate) %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.text.x = element_text(size = 10))


XJB_plots <- plot_grid(
  XJB,
  XJB_mixture_plots,
  nrow = 1,
  rel_widths = c(1, 3)  # make Mixture twice as wide
)


savePNGPDF(paste0(OUTDIR, "XJB_plots_36"), XJB_plots, 4, 10)


