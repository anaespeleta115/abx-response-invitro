# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-getFullyLostFams/100825-getFullyLostFams.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-testNullModelLostFamilies/out/"

COMPONENTS <- c("0", "1", "2", "3")
PALETTE.COMPONENTS <- c(  "cornflowerblue","coral2", "darkorange","chartreuse3")
names(PALETTE.COMPONENTS) <- COMPONENTS


curr_replicate <- 1



# --------------------------- run null model  ----------------------------


mixture_ids_36 <- unique(mixture_ASVs %>% filter(day == "036", subject1 == "XBA") %>%
                           pull(mixture))

plot_list <- list()

for (mix in mixture_ids_36){
  ids <- unlist(strsplit(mix, "\\+"))
  donor_id_long <- ids[1]
  donor_id <- str_sub(donor_id_long, 1, -5)
  recipient_id_long <- ids[2] 
  recipient_id <- str_sub(recipient_id_long, 1, -5)
  mix_pair <- paste(donor_id, recipient_id, sep = "+")
  
  
  # Get day 36 mix otus
  day36_mixture <- actual_colonizers_results %>% 
    filter(day == "036", subject1 == recipient_id, subject2 == donor_id) %>% 
    select(-colonized_day29, colonized_day36, colonized_day64)
  
  # Get day 36 recipient otus
  recipient <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id_long, day == "036", replicate == curr_replicate)
  
  # Get day 29 colonizers
  day29_colonizers <- actual_colonizers_results %>% 
    filter(day == "029", subject1 == recipient_id, subject2 == donor_id, actual_colonizer == 1)
  
  # Get day 36 colonizers
  day36_colonizers <- actual_colonizers_results %>% 
    filter(day == "036", subject1 == recipient_id, subject2 == donor_id, actual_colonizer == 1)
  
  # recipient_lost_29_36 <- recipient_lost_29_36 %>% 
  #   mutate(subject1 = str_sub(recipient_id_long, 1, -5))
  
  # Get lost otus
  recipient_lost_families <- fully_lost_families %>% 
    filter(subject1 == recipient_id)
  
  # Add a component classification for each otu
  day36_mixture_components <- day36_mixture %>%
    mutate(component = case_when(
      # Recipient takes priority over all other categories
      OTU %in% recipient$OTU ~ "1",
      
      # Colonizers only if not already in the recipient
      OTU %in% day29_colonizers$OTU & OTU %in% day36_colonizers$OTU ~ "2",
      
      # Lost strains if not in recipient or colonizers
      Family %in% recipient_lost_families$Family ~ "3",
      
      # All others
      TRUE ~ "0"
    ))
  
  
  day36_mixture_components_plot <- day36_mixture_components %>% 
    ggplot(aes(x = mixture, y = relAbundance, fill = component))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = PALETTE.COMPONENTS) +
    labs(title = "" , x = "", y = "", fill = "Component") +
    facet_wrap(~ biosample1)+
    DEFAULTS.THEME_PRINT+
    theme(legend.position = "none",
          strip.text.x = element_blank()) # remove facet labels
  
  plot_list[[mix]] <- day36_mixture_components_plot
  
}


combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 10)

savePNGPDF(paste0(OUTDIR, "day36_mix_components_fully_lost_fams"), combined_plot, 2, 9)

