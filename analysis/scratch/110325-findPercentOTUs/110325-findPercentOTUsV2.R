# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102225-gete0041colonizationProp/102225-gete0041colonizationProp.R")
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-gete0041colonization/102125-gete0041colonization.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/110325-findPercentOTUs/out/"




# ------------------------------------------ % abundance lost out of total OTUs (VERSION 1) -----------------------------------------------

e0041_control_recipients_pre_abx <- e0041_control_recipients %>% 
  filter(recipient == "pre-abx") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

e0041_control_recipients_post_abx1 <- e0041_control_recipients %>% 
  filter(recipient == "post-abx-V1") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

family_otus_pre_abx <- e0041_control_recipients_pre_abx %>% 
  group_by(Family) %>% 
  summarize(fam_otus_pre_abx = n())

family_otus_post_abx1 <- e0041_control_recipients_post_abx1 %>% 
  group_by(Family) %>% 
  summarize(fam_otus_post_abx1 = n())

family_otus <- family_otus_pre_abx %>% 
  left_join(family_otus_post_abx1, by = c("Family")) %>% 
  mutate(fam_otus_post_abx1 = ifelse(is.na(fam_otus_post_abx1), 0, fam_otus_post_abx1), 
         difference = as.numeric(fam_otus_pre_abx) - as.numeric(fam_otus_post_abx1)) %>% 
  mutate(percentage_lost = difference/fam_otus_pre_abx,
         lost = ifelse(percentage_lost > 0, 1, 0),
         gained = ifelse(percentage_lost < 0, 1, 0))

# ------------------------------------------ % of family lost (in OTU count) -----------------------------------------------



percent_lost_otus_family_plot <- family_otus %>% 
  ggplot(aes(x = fct_reorder(Family, -fam_otus_pre_abx), y = percentage_lost, fill = Family)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of family otus lost",
       fill = "Family")+
  scale_fill_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_lost_otus_fams_post_abx_v1"), percent_lost_otus_family_plot , 4, 4)



# ------------------------------------------ %  differential colonizer OTUs out of the colonizer OTUs per family (VERSION 1) -----------------------------------------------

colonizer_families <- mixture_colonization_full %>% 
  filter(colonized_post_abx_v1 == 1) %>% 
  group_by(Family) %>% 
  summarize(family_colonizer_otus = n())

diff_colonizer_families <- mixture_colonization_full %>% 
  filter(diff_colonizer_v1 == 1) %>% 
  group_by(Family) %>% 
  summarize(family_diff_colonizer_otus = n())

percentage_differential <- colonizer_families %>% 
  left_join(diff_colonizer_families, by = c("Family")) %>% 
  mutate(family_diff_colonizer_otus = ifelse(is.na(family_diff_colonizer_otus), 0, family_diff_colonizer_otus), 
         percentage = as.numeric(family_diff_colonizer_otus) / as.numeric(family_colonizer_otus) * 100,
         above_50 = ifelse(percentage > 50, 1, 0))
  

# ------------------------------------------ plot -----------------------------------------------


percent_colonizer_otus_family_plot <- percentage_differential %>% 
  ggplot(aes(x = fct_reorder(Family, -percentage), y = percentage, fill = Family, alpha = factor(above_50))) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% differential colonizers",
       fill = "Family")+
  scale_alpha_manual(values = c("0" = 0, "1" = 1))+
  scale_fill_manual(values = my_colors)+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_diff_colonizer_otus_post_abx_v1_ADJUSTED"), percent_colonizer_otus_family_plot , 3, 2)


colonizer_otus_family_plot <- colonizer_families %>% 
  ggplot(aes(x = fct_reorder(Family, -family_colonizer_otus), y = family_colonizer_otus, fill = Family)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "Number of colonizer otus",
       fill = "Family")+
  scale_fill_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "colonizer_otus_post_abx_v1"), colonizer_otus_family_plot , 4, 4)



# ------------------------------------------ %  differential colonizer OTUs out of the colonizer OTUs per family (VERSION 2) -----------------------------------------------

colonizer_families <- mixture_colonization_full %>% 
  filter(colonized_post_abx_v1 == 1) %>% 
  group_by(Family) %>% 
  summarize(family_colonizer_otus = n())

diff_colonizer_families <- mixture_colonization_full %>% 
  filter(diff_colonizer_v1 == 1) %>% 
  group_by(Family) %>% 
  summarize(family_diff_colonizer_otus = n())

percentage_differential <- colonizer_families %>% 
  left_join(diff_colonizer_families, by = c("Family")) %>% 
  mutate(family_diff_colonizer_otus = ifelse(is.na(family_diff_colonizer_otus), 0, family_diff_colonizer_otus), 
         percentage = as.numeric(family_diff_colonizer_otus) / as.numeric(family_colonizer_otus) * 100,
         above_50 = ifelse(percentage > 50, 1, 0))


# ------------------------------------------ plot -----------------------------------------------


percent_colonizer_otus_family_plot <- percentage_differential %>% 
  ggplot(aes(x = fct_reorder(Family, -percentage), y = percentage, fill = Family, alpha = factor(above_50))) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of all post-abx colonizers
       that are post-abx-only",
       fill = "Family")+
  scale_alpha_manual(values = c("0" = 0, "1" = 1))+
  scale_fill_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_diff_colonizer_otus_post_abx_v2_ADJUSTED"), percent_colonizer_otus_family_plot , 3, 2)


colonizer_otus_family_plot <- colonizer_families %>% 
  ggplot(aes(x = fct_reorder(Family, -family_colonizer_otus), y = family_colonizer_otus, fill = Family)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "Number of colonizer otus",
       fill = "Family")+
  scale_fill_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "colonizer_otus_post_abx_v2"), colonizer_otus_family_plot , 3, 3)

pre_abx_colonizer_families <- mixture_colonization_full %>% 
  filter(colonized_pre_abx == 1) %>% 
  group_by(Family) %>% 
  summarize(family_colonizer_otus = n())

pre_abx_colonizer_otus_family_plot <- pre_abx_colonizer_families %>% 
  ggplot(aes(x = fct_reorder(Family, -family_colonizer_otus), y = family_colonizer_otus, fill = Family)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "Number of colonizer otus",
       fill = "Family")+
  scale_fill_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))

savePNGPDF(paste0(OUTDIR, "pre_abx_colonizer_otus"), pre_abx_colonizer_otus_family_plot , 4, 4)


