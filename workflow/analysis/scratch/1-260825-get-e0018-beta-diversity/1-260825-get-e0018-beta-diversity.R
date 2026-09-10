source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/"


# # This script generates and plots JSD only for the autologous communities in e0018
#
# # subset phyloseq object into desired groups
# 
# # get updated metadata
# new_meta <- e0018_autologous_mixes %>%
#   select(sample, biosample_full, mix_type) %>%
#   distinct()
# 
# meta <- data.frame(sample_data(e0018_ps))
# meta$sample <- rownames(meta)
# 
# # remove old versions so updated values don't create .x/.y columns
# meta_updated <- meta %>%
#   select(-any_of(c("biosample_full", "mix_type"))) %>%
#   inner_join(new_meta, by = "sample")
# 
# # filter phyloseq object
# e0018_obj_filtered <- prune_samples(
#   meta_updated$sample,
#   e0018_ps
# )
# 
# # phyloseq expects metadata rownames = sample names
# rownames(meta_updated) <- meta_updated$sample
# meta_updated$sample <- NULL
# 
# # force identical ordering
# meta_updated <- meta_updated[
#   sample_names(e0018_obj_filtered),
#   ,
#   drop = FALSE
# ]
# 
# sample_data(e0018_obj_filtered) <- sample_data(meta_updated)
# 
# 
# 
# # Use JSD to compare same sample, different replicates, then compare different samples, different replicates.
# 
# ## UNCOMMENT THIS ONLY WHEN NECESSARY TO RE-COMPUTE JSD
# 
# distMethods <- c("jsd")
# 
# calculateBeta <- function(data, distMethod) {
#   # Calculate the distance matrix using the specified method.
#   betaRaw <- distance(data, method=distMethod)
#   
#   
#   beta <- as.matrix(betaRaw)
#   beta <- as.data.frame(beta)
#   beta$sample1 <- rownames(beta)
#   
#   beta <- beta %>%
#     pivot_longer(-sample1, names_to="sample2", values_to="value")
#   beta <- beta %>%
#     filter(sample1 != sample2) %>%
#     mutate(method=distMethod)
# }
# 
# # Calculate the distance matrix for all of the specified methods on the species abundances.
# # Combine the distance matrices for all methods.
# betaSpecies <- do.call(rbind, lapply(distMethods, function(distMethod) {
#   print(distMethod)
#   calculateBeta(e0018_obj_filtered, distMethod)
# }))
# # Export the distance matrix generated for all of the sample pairs
# # using all of the specified methods.
# write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta.txt.gz"))
# 
# 
# # Add sample IDs back in
# 
# beta_diversity <- read_delim("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/speciesBeta.txt.gz")
# 
# 
# # Join the metadata back in:
# 
# meta_updated_1 <- meta_updated %>% 
#   select(fullSample, biosample_full, mix_type, biosample1, biosample2)
# 
# beta_diversity_1 <- beta_diversity %>% 
#   left_join(meta_updated_1, by = c("sample1" = "fullSample"))
# 
# beta_diversity_1 <- beta_diversity_1 %>% 
#   rename(biosample_full_1 = biosample_full, mix_type_1 = mix_type, biosample1_1 = biosample1,
#          biosample2_1 = biosample2)
# 
# beta_diversity_final <- beta_diversity_1 %>% 
#   left_join(meta_updated_1, by = c("sample2" = "fullSample"))
# 
# 
# beta_diversity_final <- beta_diversity_final %>% 
#   rename(biosample_full_2 = biosample_full, mix_type_2 = mix_type, biosample1_2 = biosample1,
#          biosample2_2 = biosample2)
# 
# write_delim(beta_diversity_final, paste0(OUTDIR, "analogous_mixes_beta_final.txt.gz"))
# 

# ---------------------- START HERE FOR PLOTTING -----------------------

## Read final, full beta diversity dataset back in:
beta_diversity_final <- read_delim("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/analogous_mixes_beta_final.txt.gz")



# Filter out within sample comparisons???

same_mixes <- beta_diversity_final %>% 
  filter(biosample_full_1 == biosample_full_2)

diff_mixes <-  beta_diversity_final%>% 
  filter(biosample_full_1 != biosample_full_2)

# get the same mixture comparisons across replicates
diff_replicates <- same_mixes %>% 
  mutate(replicate_1 = str_sub(mix_type_1, -1), replicate_2 = str_sub(mix_type_2,-1)) %>% 
  filter(replicate_1 != replicate_2) %>% 
  mutate(combination = "same_mix_diff_replicates")

# From here, get the combination of mix_type_1 and mix_type_2 to 

relevant_comparisons <- diff_mixes %>% 
  mutate(combination = ifelse(str_sub(biosample_full_1, 2, 2) == str_sub(biosample_full_2, 2, 2), "same_household_autologous", "diff_household_autologous")) %>% 
  filter(mix_type_1 == "A-pre-1", mix_type_2 %in% c(
    "A-pre to A-post-1",
    "A-pre to A-post-2",
    "A-pre to A-post-3" # here, we've now filtered only for autologous mixture comparisons with A-pre recipients either from the same or different household
  )) %>% 
  rbind(diff_replicates %>% select(-replicate_1, -replicate_2)) %>% # remove extra columns so that they match, we already filtered for rep 1 and rep 2 to be different
  mutate(combination = fct_relevel(combination, "same_mix_diff_replicates", "same_household_autologous","diff_household_autologous"))


### ----------------- Plot boxplots across the different comparisons

p_relevant_comparisons <- relevant_comparisons %>% 
  ggplot(aes(x = combination, y = value, fill = combination)) +
  geom_boxplot()  +
  labs(x = "", y = "JSD between communities") +
  DEFAULTS.THEME_PRINT+ 
  scale_x_discrete(
    labels = c(
      "same_household_autologous" = "Same A-pre,
A-pre to A-post",
      "diff_household_autologous" = "Different A-pre,
A-pre to A-post",
      "same_mix_diff_replicates" = "Different replicates 
same mixture"
    ))+
  theme(
      legend.position = "none"
    )+
  scale_fill_paletteer_d("MetBrewer::Hiroshige")
  

p_relevant_comparisons

savePNGPDF(paste0(OUTDIR, "p_relevant_comparisons"), p_relevant_comparisons, 2.5, 3)



### ----------------- Plot heatmap to compare autologous mixtures and their avg. JSD across replicates

p_autologous_heatmap <- relevant_comparisons |> 
  filter(combination == "same_household_autologous") |> 
  group_by(biosample_full_1, biosample_full_2) |> 
  summarize(avg_value = mean(value)) |> 
  ggplot(aes(x = biosample_full_2, y = biosample_full_1, fill = avg_value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A",  limits = c(0, 1)) +
  labs(x = "A-pre to A-post", y = "A-pre", fill = "Avg. JSD",
       title = "")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  DEFAULTS.THEME_PRINT

p_autologous_heatmap 

savePNGPDF(paste0(OUTDIR, "p_autologous_heatmap"), p_autologous_heatmap, 2.5, 3)



