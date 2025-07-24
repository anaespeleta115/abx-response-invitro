source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getRecipient-MixtureRichness/072125-getRecipient-MixtureRichness.R")
source("C:/abx-response-invitro/analysis/scratch/072225-getColonizationProportion/072225-getColonizationProportion.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072125-compareRichness-Colonization/out/"

# We need the colonization_results and the recipient_richness datasets. From there, we have to join them by biosample1 

total_colonization <- actual_colonizers_results %>% 
  group_by(mixture) %>% 
  summarize(total_colonizers = sum(actual_colonizer)) %>% 
  mutate(recipient = str_sub(mixture, -7, -1),
         donor = str_sub(mixture, 1, -9))


combined_richness <- left_join(total_colonization, recipient_richness,
                               by = c("recipient" = "sample") # match recipient to sample
) %>% 
  select(-community)  %>% 
  mutate(recipient = str_sub(recipient, -7, -5),
         donor = str_sub(mixture, 1, -9))


p_colonization_richness <- combined_richness %>% 
  ggplot(aes(x = species_richness, y = total_colonizers, color = factor(day))) + 
  geom_point()+
  geom_text_repel(aes(label = recipient), size = 2.5, max.overlaps = 20)+
  labs(
    x = "Species Richness of Recipient",
    y = "Number of Colonizers in Mixture",
    color = "Study Day"
  )+
  facet_wrap(~donor)+
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "colonizationProp-richness"), p_colonization_richness, 4, 8)


