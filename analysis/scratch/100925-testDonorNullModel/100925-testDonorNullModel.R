# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-getFullyLostFams/100825-getFullyLostFams.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100925-testDonorNullModel/out/"

COMPONENTS <- c("0", "1", "2", "3", "4", "5")
PALETTE.COMPONENTS <- c(  "darkgrey","coral2", "darkorange","chartreuse3", "#ffc425", "darkorchid")
names(PALETTE.COMPONENTS) <- COMPONENTS


curr_replicate <- 1


# --------------------------- Try out one community (XBB + XBA) ----------------------------

XBB_donor_day36 <- single_donor_ASVs %>% 
  filter(biosample1 == "XBB-029")

# Get day 29 colonizers
XBAXBB_day29_colonizers <- actual_colonizers_results %>% 
  filter(day == "029", subject2 == "XBB", actual_colonizer == 1)

# Get day 36 colonizers
XBAXBB_day36_colonizers <- actual_colonizers_results %>% 
  filter(day == "036", subject2 == "XBB", actual_colonizer == 1)

# Get day 36 recipient otus
XBA_recipient <- recipient_ASVs %>% 
  filter(biosample1 == "XBA-036")

# Get lost otus
XBA_lost_strains <- fully_lost_families %>% 
  filter(subject1 == "XBA")

# Add a component classification for each otu
XBAXBB_donor_components <- XBB_donor_day36 %>% 
  mutate(component = case_when(
    OTU %in% XBA_recipient$OTU ~ "1", # OTUs present in post-abx recipient
    OTU %in% XBAXBB_day29_colonizers$OTU ~ "2", # OTUs present as colonizers in day 29 AND day 36 mixtures
    Family %in% XBA_lost_strains$Family ~ "3", # OTUs lost post-abx in the recipient
    TRUE ~ "0" # "other"
  ))


# --------------------------- PLOT XBB WITH COMPONENTS  ----------------------------

XBAXBB_day36_mixture_components_plot <- XBAXBB_donor_components %>% 
  ggplot(aes(x = biosample1, y = relAbundance, fill = component))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = PALETTE.COMPONENTS) +
  labs(title = "" , x = "Donor", y = "Relative abundance", fill = "Component") +
  facet_wrap(~ biosample1)+
  DEFAULTS.THEME_PRINT



# --------------------------- run model on all mixtures ----------------------------


mixture_ids_36 <- unique(mixture_ASVs %>% filter(day == "036") %>%
                           pull(mixture))

all_donors_components <- tibble()

plot_list <- list()

for (mix in mixture_ids_36){
  ids <- unlist(strsplit(mix, "\\+"))
  donor_id_long <- ids[1]
  donor_id <- str_sub(donor_id_long, 1, -5)
  recipient_id_long <- ids[2] 
  recipient_id <- str_sub(recipient_id_long, 1, -5)
  mix_pair <- paste(donor_id, recipient_id, sep = "+")
  
  # Get donor otus
  donor_day36 <- single_donor_ASVs %>% 
    filter(biosample1 == donor_id_long)
  
  # # Get day 36 mix otus
  # day36_mixture <- actual_colonizers_results %>% 
  #   filter(day == "036", subject1 == recipient_id, subject2 == donor_id) %>% 
  #   select(-colonized_day29, colonized_day36, colonized_day64)
  
  # Get day 36 recipient otus
  recipient <- recipient_ASVs %>% 
    filter(biosample1 == recipient_id_long, day == "036", replicate == curr_replicate)
  
  # Get day 29 colonizers
  day29_colonizers <- actual_colonizers_results %>% 
    filter(day == "029", subject1 == recipient_id, subject2 == donor_id, actual_colonizer == 1)
  
  # Get day 36 colonizers
  day36_colonizers <- actual_colonizers_results %>% 
    filter(day == "036", subject1 == recipient_id, subject2 == donor_id, actual_colonizer == 1)
  
  # Get lost otus
  recipient_lost_families <- fully_lost_families %>% 
    filter(subject1 == recipient_id)
  
  # Add a component classification for each otu
  day36_donor_components <- donor_day36 %>%
    mutate(component = case_when(
      # Recipient takes priority over all other categories
      OTU %in% recipient$OTU ~ "1",
      
      # Colonizers only if not already in the recipient
      !(OTU %in% recipient$OTU) & OTU %in% day29_colonizers$OTU & OTU %in% day36_colonizers$OTU ~ "2",
      
      # Lost strains if not in recipient but in day 36 colonizers
      !(OTU %in% recipient$OTU) & Family %in% recipient_lost_families$Family & OTU %in% day36_colonizers$OTU ~ "3",
      
      # Lost strains in donor not in day 36 colonizers
      !(OTU %in% recipient$OTU) & Family %in% recipient_lost_families$Family & !(OTU %in% day36_colonizers$OTU) ~ "4",
      
      # Non-lost strains in donor that did colonize day 36 (in theory these should be minimal)
      !(OTU %in% recipient$OTU) & !(Family %in% recipient_lost_families$Family) & OTU %in% day36_colonizers$OTU ~ "5",
      
      # All others: non-lost strains in donor that did not colonize (these should be more than the purple)
      TRUE ~ "0"
    ))
  
  day36_donor_components <- day36_donor_components %>% 
    mutate(mixture = mix)
  
  day36_donor_components_plot <- day36_donor_components %>% 
    ggplot(aes(x = biosample1, y = relAbundance, fill = component))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = PALETTE.COMPONENTS) +
    labs(title = "" , x = "", y = "", fill = "Component") +
    facet_wrap(~ biosample1)+
    DEFAULTS.THEME_PRINT+
    theme(legend.position = "none",
          strip.text.x = element_blank()) # remove facet labels
  
  plot_list[[mix]] <- day36_donor_components_plot
  
  all_donors_components <- bind_rows(all_donors_components, day36_donor_components)
  
}


combined_plot <- wrap_plots(plot_list, nrow = 4, ncol = 10)

savePNGPDF(paste0(OUTDIR, "day36_mix_components_fully_lost_fams_donors_XKA"), combined_plot, 6, 9)

