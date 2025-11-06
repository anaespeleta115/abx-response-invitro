# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/110325-findPercentOTUs/out/"


# This version is not the most updated. 

# ------------------------------------------ % abundance lost out of total OTUs (VERSION 1) -----------------------------------------------

e0041_control_recipients_pre_abx <- e0041_control_recipients %>% 
  filter(recipient == "pre-abx") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

e0041_control_recipients_post_abx1 <- e0041_control_recipients %>% 
  filter(recipient == "post-abx-V1") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

family_abundance_pre_abx <- e0041_control_recipients_pre_abx %>% 
  group_by(Family) %>% 
  summarize(fam_abundance_pre_abx = sum(relAbundance))

family_abundance_post_abx1 <- e0041_control_recipients_post_abx1 %>% 
  group_by(Family) %>% 
  summarize(fam_abundance_post_abx1 = sum(relAbundance))


family_abundance <- family_abundance_pre_abx %>% 
  left_join(family_abundance_post_abx1, by = c("Family")) %>% 
  mutate(fam_abundance_post_abx1 = ifelse(is.na(fam_abundance_post_abx1), 0, fam_abundance_post_abx1), 
         difference = as.numeric(fam_abundance_pre_abx) - as.numeric(fam_abundance_post_abx1)) %>% 
  mutate(percentage_lost = difference/fam_abundance_pre_abx,
         lost = ifelse(percentage_lost > 0, 1, 0),
         gained = ifelse(percentage_lost < 0, 1, 0))

# ------------------------------------------ % of family lost (in OTU count) -----------------------------------------------

lost_family_abundance <- family_abundance %>% 
  filter(lost == 1) 

percent_lost_family_plot <- lost_family_abundance %>% 
  ggplot(aes(x = fct_reorder(Family, -fam_abundance_pre_abx), y = percentage_lost, fill = fam_abundance_pre_abx)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of family abundance lost",
       fill = "Pre-abx family abundance")+
  scale_fill_viridis_c(option = "magma")+ # limits = c(0.0001, 0.3)
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_lost_fams_post_abx_v1"), percent_lost_family_plot , 4, 4)



# ------------------------------------------ % family differential colonizers out of total family abundance (VERSION 1) -----------------------------------------------

e0041_control_recipients_pre_abx <- e0041_control_recipients %>% 
  filter(recipient == "pre-abx") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

e0041_control_recipients_post_abx1 <- e0041_control_recipients %>% 
  filter(recipient == "post-abx-V1") %>% 
  select(recipient, Family, OTU, relAbundance) %>% 
  distinct()

family_abundance_pre_abx <- e0041_control_recipients_pre_abx %>% 
  group_by(Family) %>% 
  summarize(fam_abundance_pre_abx = sum(relAbundance))

family_abundance_post_abx1 <- e0041_control_recipients_post_abx1 %>% 
  group_by(Family) %>% 
  summarize(fam_abundance_post_abx1 = sum(relAbundance))


family_abundance <- family_abundance_pre_abx %>% 
  left_join(family_abundance_post_abx1, by = c("Family")) %>% 
  mutate(fam_abundance_post_abx1 = ifelse(is.na(fam_abundance_post_abx1), 0, fam_abundance_post_abx1), 
         difference = as.numeric(fam_abundance_pre_abx) - as.numeric(fam_abundance_post_abx1)) %>% 
  mutate(percentage_lost = difference/fam_abundance_pre_abx,
         lost = ifelse(percentage_lost > 0, 1, 0),
         gained = ifelse(percentage_lost < 0, 1, 0))

# ------------------------------------------ % of family lost (in abundance) -----------------------------------------------

lost_family_abundance <- family_abundance %>% 
  filter(lost == 1) 

percent_lost_family_plot <- lost_family_abundance %>% 
  ggplot(aes(x = fct_reorder(Family, -fam_abundance_pre_abx), y = percentage_lost, fill = fam_abundance_pre_abx)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "% of family abundance lost",
       fill = "Pre-abx family abundance")+
  scale_fill_viridis_c(option = "magma")+ # limits = c(0.0001, 0.3)
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 6, angle = 90))


savePNGPDF(paste0(OUTDIR, "percentage_lost_fams_post_abx_v1"), percent_lost_family_plot , 4, 4)


# FOR DIFFERENTIAL COLONIZATION: for every colonizer family in the mixture community, what percent of their OTUs were differential
# colonizers??


