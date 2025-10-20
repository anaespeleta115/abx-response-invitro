# Load data

source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/100925-testDonorNullModel/100925-testDonorNullModel.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100925-testDonorNullModel/out/"




# --------------------------- CALCULATE COMPONENT 1 PERCENTAGE  ----------------------------

donor_components_percentages <- all_donors_components  %>% 
  group_by(subject1, subject2, component) %>% 
  summarize(total_otus_component = sum(relAbundance)) %>% 
  pivot_wider(names_from = component, values_from = total_otus_component) %>% 
  rename(component_0 = "0", component_1 = "1", component_2 = "2", component_3 = "3") %>% 
  replace_na(list(component_0 = 0, component_3 = 0)) %>% 
  mutate(non_recipient = sum(component_0, component_2, component_3), recipient = sum(component_1), 
         total = sum(component_0, component_1, component_2, component_3), null_model = sum(component_1, component_2, component_3))


component_1_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_1_percentage = recipient/total * 100) 

component_1_percentage_plot <- ggplot(component_1_percentage, aes(x = fct_reorder(subject1, -component_1_percentage), y = component_1_percentage)) +
  geom_boxplot(fill = "#ffc425") +
  # geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Recipient's mixtures", y = "Percentage of colonizers belonging to component 1") +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "component_1_percentage_plot_v3_recipient"), component_1_percentage_plot, 3, 3)


# Identify unexpected colonizers (component 5)

all_donors_components <- all_donors_components %>% mutate(subject1 = str_sub(mixture, 9, -1))

unexpected_colonizers <- all_donors_components %>% 
  filter(component == 5) %>% 
  group_by(subject1, Family) %>% 
  summarize(fam_counts = n())

lost_taxa_colonizers <- all_donors_components %>% 
  filter(component == 3) %>% 
  group_by(mixture, Family) %>% 
  summarize(fam_counts = n())

uninvolved_taxa <- all_donors_components %>% 
  filter(component == 0) %>% 
  group_by(subject1, Family) %>% 
  summarize(fam_counts = n())
  