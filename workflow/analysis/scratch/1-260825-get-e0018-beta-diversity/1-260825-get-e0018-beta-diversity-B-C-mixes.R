source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/"


# # This script generates and plots JSD for all communities in e0018
# 
# # START HERE IF NEEDING TO UPDATE METADATA
# # subset phyloseq object into desired groups
# 
# # get all relevant metadata
# new_meta <- e0018_data_full %>%
#   filter(!is.na(biosample_full)) %>%
#   select(sample, biosample_full, mixture) %>%
#   distinct()
# 
# # extract metadata from phyloseq object
# meta <- data.frame(sample_data(e0018_ps))
# meta$sample <- rownames(meta)
# 
# # meta updated now contains metadata from e0018_data_full
# meta_updated <- meta %>%
#   select(-any_of(c("biosample_full"))) %>%
#   inner_join(new_meta, by = "sample")
# 
# # use meta_updated to filter phyloseq object
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
# write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta_full.txt.gz"))
# 
# 
# # Add sample IDs back in
# 
# beta_diversity <- read_delim("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/speciesBeta_full.txt.gz")
# 
# 
# # Join the metadata back in:
# 
# meta_updated_1 <- meta_updated %>%
#   select(fullSample, biosample_full, biosample1, biosample2, mixture, inoculationReplicate, ratio)
# 
# beta_diversity_1 <- beta_diversity %>%
#   left_join(meta_updated_1, by = c("sample1" = "fullSample"))
# 
# beta_diversity_1 <- beta_diversity_1 %>%
#   rename(biosample_full_1 = biosample_full, biosample1_1 = biosample1,
#          biosample2_1 = biosample2, mixture_1 = mixture, inoculationReplicate_1 = inoculationReplicate, ratio_1 = ratio)
# 
# beta_diversity_final <- beta_diversity_1 %>%
#   left_join(meta_updated_1, by = c("sample2" = "fullSample"))
# 
# 
# beta_diversity_final <- beta_diversity_final %>%
#   rename(biosample_full_2 = biosample_full, biosample1_2 = biosample1,
#          biosample2_2 = biosample2, mixture_2 = mixture, inoculationReplicate_2 = inoculationReplicate, ratio_2 = ratio)
# 
# write_delim(beta_diversity_final, paste0(OUTDIR, "all_mixes_beta_final.txt.gz"))
# 

# ---------------------- START HERE FOR PLOTTING !!! -----------------------

## Read final, full beta diversity dataset back in:
beta_diversity_final <- read_delim("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/out/all_mixes_beta_final.txt.gz") %>%
  filter(!is.na(biosample_full_2), !is.na(biosample_full_1))


# add a combination column to represent the comparison and filter for only the relevant comparisons
beta_diversity_final_combinations <- beta_diversity_final %>% 
  mutate(combination = paste0(mixture_1,"+", mixture_2)) %>% 
  filter(combination %in% c("A-pre+A-post", "A-pre+A-pre to A-post", "A-pre+B to A-post", "A-pre+C to A-post", "A-pre+A-post")) %>% 
  mutate(full_combination = paste0(biosample_full_1,"+",biosample_full_2))


# add back in the same-mix-replicate comparison

same_mixes <- beta_diversity_final %>% 
  filter(biosample_full_1 == biosample_full_2)

diff_replicates <- same_mixes %>%
  filter(inoculationReplicate_1 != inoculationReplicate_2, ratio_2 == "1:100", ratio_1 == "1:100") %>%
  mutate(combination = "same_mix_diff_replicates", full_combination = paste0(biosample_full_1,"+",biosample_full_2))
  


# From here, get the combination of the full biosamples and divide them into same recipient, diff recipient or same community.

relevant_comparisons <- beta_diversity_final_combinations %>%
  mutate(type = ifelse(str_sub(biosample1_1, 2, 2) == str_sub(biosample2_2, 2, 2), "same_recipient", "diff_recipient")) %>%
  filter(ratio_2 == "1:100" | (is.na(ratio_1) & is.na(ratio_2))) %>% 
  mutate(type = ifelse(combination == "A-pre+A-post" & str_sub(biosample1_1, 2, 2) == str_sub(biosample1_2, 2, 2), "same_recipient", type)) %>% 
  bind_rows(diff_replicates %>% select(-inoculationReplicate_1, -inoculationReplicate_2) %>% mutate(type = "same_community")) %>% 
  mutate(combination = fct_relevel(combination, "same_mix_diff_replicates", "A-pre+A-pre to A-post", "A-pre+B to A-post", "A-pre+C to A-post", "A-pre+A-post"))



p_relevant_comparisons_full <- relevant_comparisons %>% 
  filter(type == "same_recipient" | type == "same_community") %>% 
  ggplot(aes(x = combination, y = value, fill = combination)) +
  geom_boxplot()  +
  labs(x = "", y = "JSD between communities") +
  DEFAULTS.THEME_PRINT+
  scale_x_discrete(
    labels = c(
      "same_mix_diff_replicates" = "Replicates",
      "A-pre+A-pre to A-post" = "A-pre vs.
A-pre to A-post",
      "A-pre+A-post" = "A-pre vs.
      A-post",
      "A-pre+B to A-post" = "A-pre vs.
B to A-post",
      "A-pre+C to A-post" = "A-pre vs.
C to A-post"
    ))+
  # missing A-pre vs. A-post and replicates
  theme(
      legend.position = "none"
    )+
  scale_fill_paletteer_d("MetBrewer::Hiroshige")
  

p_relevant_comparisons_full

savePNGPDF(paste0(OUTDIR, "p_relevant_comparisons_full"), p_relevant_comparisons_full, 2.5, 3)


# --------------- plot heatmap for A-pre+A-post comparisons

recipient_comparisons <- relevant_comparisons |> 
  filter(combination == "A-pre+A-post", type == "same_recipient")
  


p_recipient_comparisons <- recipient_comparisons |> 
  group_by(biosample_full_1, biosample_full_2) |> 
  ggplot(aes(x = biosample_full_2, y = biosample_full_1, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako",  limits = c(0, 1)) +
  labs(x = "A-pre", y = "A-post", fill = "JSD",
       title = "")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  DEFAULTS.THEME_PRINT

p_recipient_comparisons

savePNGPDF(paste0(OUTDIR, "p_recipient_comparisons"), p_recipient_comparisons, 2.5, 3)










