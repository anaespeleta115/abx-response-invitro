# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/103125-plotFamRecovery/out/"



# Get fam relAbundance pre-abx

fam_abundance_pre_abx <- e0041_control_recipients %>% 
  filter(recipient == "pre-abx") %>% 
  group_by(Family) %>% 
  summarize(pre_abx_fam_relAbundance = sum(relAbundance))


# Get fam relAbundance post-abx V1

fam_abundance_post_abx_V1 <- e0041.A %>% 
  filter(recipient == "post-abx-V1") %>% 
  group_by(mixture, Family) %>% 
  summarize(post_abx1_fam_relAbundance = sum(relAbundance)) 

# Get fam relAbundance post-abx V2

fam_abundance_post_abx_V2 <- e0041.A %>% 
  filter(recipient == "post-abx-V2") %>% 
  group_by(mixture, Family) %>% 
  summarize(post_abx2_fam_relAbundance = sum(relAbundance)) 

# Combine data

combined1 <- fam_abundance_pre_abx %>% 
  full_join(fam_abundance_post_abx_V1, by="Family") %>% 
  mutate(donor = str_sub(mixture, 13, -1)) %>% 
  group_by(Family, pre_abx_fam_relAbundance) %>% 
  summarize(mean_post_abx1_fam_relAbundance = mean(post_abx1_fam_relAbundance)) %>% 
  replace_na(list(pre_abx_fam_relAbundance = 0, mean_post_abx1_fam_relAbundance = 0))

# Plot recovery V1

recovery_plot1 <- combined1 %>% 
  # filter(recipient_fam_relAbundance_day_029 != 1e-04, recipient_fam_relAbundance_day_036 == 1e-04) %>%
  ggplot(aes(x = log10(pre_abx_fam_relAbundance), y = log(mean_post_abx1_fam_relAbundance) , color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  labs(x = "")+
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "recovery_plot1"), recovery_plot1, 3, 4)

# Plot recovery V2 :  TODO: THE SAME AS ABOVE BUT FOR V2

combined2 <- fam_abundance_pre_abx %>% 
  full_join(fam_abundance_post_abx_V2, by="Family") %>% 
  mutate(donor = str_sub(mixture, 13, -1)) %>% 
  group_by(Family, pre_abx_fam_relAbundance) %>% 
  summarize(mean_post_abx2_fam_relAbundance = mean(post_abx2_fam_relAbundance)) %>% 
  replace_na(list(pre_abx_fam_relAbundance = 0, mean_post_abx2_fam_relAbundance = 0))


recovery_plot2 <- combined2 %>% 
  ggplot(aes(x = log10(pre_abx_fam_relAbundance), y = log10(mean_post_abx2_fam_relAbundance) , color = Family))+
  geom_point(size = 2, alpha = 0.5)+
  labs(x = "")+
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")+
  DEFAULTS.THEME_PRINT

savePNGPDF(paste0(OUTDIR, "recovery_plot2"), recovery_plot2, 3, 4)


