# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-testNullModelLostFamilies/100825-testNullModelLostFamilies.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100725-calculateComponentPercentages/out/model3/"

curr_replicate <- 1



# --------------------------- CALCULATE COMPONENT 1 PERCENTAGE  ----------------------------

day36_mixture_components_percentages <- all_mix_day36_components  %>% 
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


# --------------------------- CALCULATE PERCENTAGE OF COLONIZERS THAT BELONG TO THE LOST TAXA FAMILIES  ----------------------------

component_3_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_3_percentage = component_3/non_recipient * 100) 

component_3_percentage_plot <- ggplot(component_3_percentage, aes(x = fct_reorder(subject1, -component_3_percentage), y = component_3_percentage)) +
  geom_boxplot(fill = "#ffc425") +
  # geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Recipient's mixtures", y = "Percentage of colonizers belonging to component 3") +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "component_3_percentage_plot_v3_recipient"), component_3_percentage_plot, 4, 3)

# --------------------------- CALCULATE PERCENTAGE OF COMPONENT 0  ----------------------------


component_0_percent <- day36_mixture_components_percentages %>% 
  ungroup() %>% 
  mutate(component_0_percentage = null_model/total * 100) %>% 
  summarize(avg = mean(component_0_percentage))


component_0_percentage_plot <- ggplot(component_0_percent, aes(x = fct_reorder(subject2, -component_0_percentage), y = component_0_percentage)) +
  # geom_boxplot(fill = "#ffc425") +
  geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Donor's mixtures", y = "Percentage of OTUs explained by model") +
  DEFAULTS.THEME_PRINT
# theme(axis.text.x  = element_blank(),
#       axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "component_0_percentage_plot_v3_donor"), component_0_percentage_plot, 3, 3)


# --------------------------- CALCULATE ONLY PERCENTAGE 3  ----------------------------

component_3_percentage <- day36_mixture_components_percentages %>% 
  mutate(component_3_percentage = component_3/total * 100) 

component_3_percentage_plot <- ggplot(component_3_percentage, aes(x = fct_reorder(subject2, -component_3_percentage), y = component_3_percentage)) +
  # geom_boxplot(fill = "#ffc425") +
  geom_boxplot(fill = "#d11141") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Donor's mixtures", y = "Percentage of OTUs belonging to component 3") +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "component_3_percentage_plot_ONLY3_donor"), component_3_percentage_plot, 4, 3)


