source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getRecipient-MixtureRichness/072125-getRecipient-MixtureRichness.R")
source("C:/abx-response-invitro/analysis/scratch/072125-getColonization/072125-getColonization.R")
source("C:/abx-response-invitro/analysis/scratch/072225-getColonizationProportion/072225-getColonizationProportion.R")



# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/072125-compareRichness-Colonization/out/"

# We need the colonization_results and the recipient_richness datasets. From there, we have to join them by biosample1 

# total_colonization <- colonization_prop_results %>% 
#   group_by(mixture) %>% 
#   summarize(total_colonizers = sum(actual_colonizer)) %>% 
#   mutate(recipient = str_sub(mixture, -7, -1),
#          donor = str_sub(mixture, 1, -9))


combined_richness <- left_join(colonization_prop_results, recipient_richness %>% select(sample, species_richness),
                               by = c("biosample1" = "sample") # match recipient to sample
) %>% filter(donor == "XHB")


p_colonization_richness <- combined_richness %>%
  ggplot(aes(x = species_richness, y = prop_colonizers, color = factor(day))) + 
  geom_point(aes(shape = recipient), size = 2)+
  # geom_text_repel(size = 2, max.overlaps = 20)+
  labs(
    x = "Recipient species richness",
    y = "Colonization proportion",
    color = "Study day",
    shape = "Recipient"
  )+
  # facet_wrap(~donor)+
  scale_color_manual(values = PALETTE.DAY)+
  DEFAULTS.THEME_PRINT +
  theme(legend.key.height = unit(7, "pt"),
        legend.key.spacing = unit(2, "pt"))


savePNGPDF(paste0(OUTDIR, "colonizationProp-richness"), p_colonization_richness, 1.5, 2.5)




# Maybe try doing the same with donor richness??


combined_richness_donor <- left_join(colonization_prop_results, donor_richness %>% select(sample, species_richness),
                               by = c("biosample2" = "sample") # match recipient to sample
)


p_colonization_richness_donor <- combined_richness_donor %>% filter(recipient == "XBA") %>%
  ggplot(aes(x = species_richness, y = prop_colonizers, color = factor(day))) + 
  geom_point(size = 2)+
  geom_text_repel(aes(label = donor), size = 2, max.overlaps = 20)+
  labs(
    x = "Donor species richness",
    y = "Colonization proportion",
    color = "Study day"
  )+
  # facet_wrap(~recipient)+
  scale_color_manual(values = PALETTE.DAY)+
  DEFAULTS.THEME_PRINT +
  theme(legend.key.height = unit(7, "pt"),
        legend.key.spacing = unit(2, "pt"))


savePNGPDF(paste0(OUTDIR, "colonizationProp-richness-donor"), p_colonization_richness_donor, 2, 4)
