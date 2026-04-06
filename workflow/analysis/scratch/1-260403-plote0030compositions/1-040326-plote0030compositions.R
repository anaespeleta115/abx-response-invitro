source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260403-loade0030data/1-260403-loade0030data.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260403-plote0030compositions/out/"


# ------------- Plot XBA-029 mixtures


# Filter to include only one subject
e0030_XBA <- e0030 %>% 
  filter(biosample1 == "XBA-029", replicate == curr_replicate) %>% 
  mutate(mixture = case_when(household2 == "blank" & household3 == "blank" & well == "A11" ~ mixture,
                             household2 == "blank" & household3 == "blank" & well == "C11" ~ paste0(mixture, ".2"),
                             # household2 == "blank" & household3 == "blank" & well == "E11" ~ paste0(mixture, ".3"),
                             TRUE ~ mixture)) %>%
  filter(biosample2 != "B-mix", biosample3 != "B-mix") 

# Add a group label:
e0030_XBA <- e0030_XBA %>%
  mutate(group = case_when(
    str_detect(household2, "\\+") ~ household2,
    TRUE ~ paste0(household2, "+", household3)
  )) %>%
  arrange(group, mixture) %>%
  mutate(mixture = factor(mixture, levels = unique(mixture)))

# Plot compositions with recipient composition
XBA_mixes_compositions <- e0030_XBA %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

XBA_mixes_compositions


savePNGPDF(paste0(OUTDIR, "XBA-029_e0030_compositions"), XBA_mixes_compositions, 5, 5)



# ------------- Plot XBA-036 mixtures


# Filter to include only one subject
e0030_XBA_036 <- e0030 %>% 
  filter(biosample1 == "XBA-036", replicate == curr_replicate) %>% 
  mutate(mixture = case_when(household2 == "blank" & household3 == "blank" & well == "A11" ~ mixture,
                             household2 == "blank" & household3 == "blank" & well == "C11" ~ paste0(mixture, ".2"),
                             household2 == "blank" & household3 == "blank" & well == "E11" ~ paste0(mixture, ".3"),
                             TRUE ~ mixture)) %>%
  filter(biosample2 != "B-mix", biosample3 != "B-mix") 

# Add a group label:
e0030_XBA_036 <- e0030_XBA_036 %>%
  mutate(group = case_when(
    str_detect(household2, "\\+") ~ household2,
    TRUE ~ paste0(household2, "+", household3)
  )) %>%
  arrange(group, mixture) %>%
  mutate(mixture = factor(mixture, levels = unique(mixture)))

# Plot compositions with recipient composition
XBA_mixes_compositions_036 <- e0030_XBA_036 %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

XBA_mixes_compositions_036


savePNGPDF(paste0(OUTDIR, "XBA-036_e0030_compositions"), XBA_mixes_compositions_036, 5, 5)



# ------------- Plot XEA-029 mixtures


# Filter to include only one subject
e0030_XEA <- e0030 %>% 
  filter(biosample1 == "XEA-029", replicate == curr_replicate) %>% 
  mutate(mixture = case_when(household2 == "blank" & household3 == "blank" & well == "A11" ~ mixture,
                             household2 == "blank" & household3 == "blank" & well == "C11" ~ paste0(mixture, ".2"),
                             # household2 == "blank" & household3 == "blank" & well == "E11" ~ paste0(mixture, ".3"),
                             TRUE ~ mixture)) %>%
  filter(biosample2 != "B-mix", biosample3 != "B-mix") 

# Add a group label:
e0030_XEA <- e0030_XEA %>%
  mutate(group = case_when(
    str_detect(household2, "\\+") ~ household2,
    TRUE ~ paste0(household2, "+", household3)
  )) %>%
  arrange(group, mixture) %>%
  mutate(mixture = factor(mixture, levels = unique(mixture)))

# Plot compositions with recipient composition
XEA_mixes_compositions <- e0030_XEA %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

XEA_mixes_compositions


savePNGPDF(paste0(OUTDIR, "XEA-029_e0030_compositions"), XEA_mixes_compositions, 5, 5)




# ------------- Plot XEA-036 mixtures


# Filter to include only one subject
e0030_XEA_036 <- e0030 %>% 
  filter(biosample1 == "XEA-036", replicate == curr_replicate) %>% 
  mutate(mixture = case_when(household2 == "blank" & household3 == "blank" & well == "A11" ~ mixture,
                             household2 == "blank" & household3 == "blank" & well == "G11" ~ paste0(mixture, ".2"),
                             household2 == "blank" & household3 == "blank" & well == "E11" ~ paste0(mixture, ".3"),
                             TRUE ~ mixture)) %>%
  filter(biosample2 != "B-mix", biosample3 != "B-mix") 

# Add a group label:
e0030_XEA_036 <- e0030_XEA_036 %>%
  mutate(group = case_when(
    str_detect(household2, "\\+") ~ household2,
    TRUE ~ paste0(household2, "+", household3)
  )) %>%
  arrange(group, mixture) %>%
  mutate(mixture = factor(mixture, levels = unique(mixture)))

# Plot compositions with recipient composition
XEA_mixes_compositions_036 <- e0030_XEA_036 %>% 
  ggplot(aes(x = mixture, y = relAbundance, fill = Family))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )+
  DEFAULTS.THEME_PRINT

XEA_mixes_compositions_036


savePNGPDF(paste0(OUTDIR, "XEA-036_e0030_compositions"), XEA_mixes_compositions_036, 5, 5)


