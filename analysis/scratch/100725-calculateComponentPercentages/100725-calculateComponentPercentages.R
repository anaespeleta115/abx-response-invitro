# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100725-calculateComponentPercentages/out/"

COMPONENTS <- c("0", "1", "2", "3")
PALETTE.COMPONENTS <- c(  "cornflowerblue","coral2", "darkorange","chartreuse3")
names(PALETTE.COMPONENTS) <- COMPONENTS


curr_replicate <- 1


# --------------------------- DAY 36 MIX WITH COMPONENTS  ----------------------------


mixture_ids_36 <- unique(mixture_ASVs %>% filter(day == "036") %>%
                           pull(mixture))

# Initialize final dataset
all_mix_day36_components <- tibble()


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
  
  recipient_lost_29_36 <- recipient_lost_29_36 %>% 
    mutate(subject1 = str_sub(recipient_id_long, 1, -5))
  
  # Get lost otus
  lost_strains <- recipient_lost_29_36 %>% 
    filter(subject1 == recipient_id) # fix this!!!
  
  # Add a component classification for each otu
  day36_mixture_components <- day36_mixture %>% 
    mutate(component = case_when(
      OTU %in% recipient$OTU ~ "1", # OTUs present in post-abx recipient
      OTU %in% day29_colonizers$OTU & OTU %in% day36_colonizers$OTU ~ "2", # OTUs present as colonizers in day 29 AND day 36 mixtures
      Family %in% lost_strains$Family ~ "3", # Families lost post-abx in the recipient
      TRUE ~ "0" # "other"
    ))
  
  all_mix_day36_components <- bind_rows(all_mix_day36_components, day36_mixture_components)

}

# --------------------------- CALCULATE COMPONENT 1 PERCENTAGE  ----------------------------

day36_mixture_components_percentages <- all_mix_day36_components  %>% 
  group_by(subject1, subject2, component) %>% 
  summarize(total_otus_component = sum(relAbundance)) %>% 
  pivot_wider(names_from = component, values_from = total_otus_component) %>% 
  rename(component_0 = "0", component_1 = "1", component_2 = "2", component_3 = "3") %>% 
  replace_na(list(component_0 = 0)) %>% 
  mutate(non_recipient = sum(component_0, component_2, component_3), recipient = sum(component_1), 
         total = sum(component_0, component_1, component_2, component_3), null_model = sum(component_1, component_2, component_3))


component_1_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_1_percentage = recipient/total * 100) 

component_1_percentage_plot <- ggplot(component_1_percentage, aes(x = subject1, y = component_1_percentage)) +
  geom_boxplot(fill = "#ffc425") +
  # geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Recipient's mixtures", y = "Percentage of colonizers belonging to component 1") +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "component_1_percentage_plot"), component_1_percentage_plot, 3, 3)

# --------------------------- CALCULATE PERCENTAGE OF COLONIZERS THAT BELONG TO THE LOST TAXA FAMILIES  ----------------------------

component_3_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_3_percentage = component_3/non_recipient * 100) 

component_3_percentage_plot <- ggplot(component_3_percentage, aes(x = subject1, y = component_3_percentage)) +
  geom_boxplot(fill = "#ffc425") +
  # geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Recipient's mixtures", y = "Percentage of colonizers belonging to component 3") +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "component_3_percentage_plot"), component_3_percentage_plot, 4, 3)

# --------------------------- CALCULATE PERCENTAGE OF COMPONENT 0  ----------------------------


component_0_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_0_percentage = null_model/total * 100)
  # group_by(subject2) %>% 
  # summarize(avg_per_donor = mean(component_0_percentage))


component_0_percentage_plot <- ggplot(component_0_percentage, aes(x = subject1, y = component_0_percentage)) +
  geom_boxplot(fill = "#ffc425") +
  # geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(80, 100)) +
  labs(x = "Recipient's mixtures", y = "Percentage of OTUs belonging to component 0") +
  DEFAULTS.THEME_PRINT
  # theme(axis.text.x  = element_blank(),
  #       axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "component_0_percentage_plot"), component_0_percentage_plot, 3, 3)







