# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/100825-getFullyLostFams/out/"

curr_replicate <- 1

recipient_ASVs <- recipient_ASVs %>%
  filter(replicate == curr_replicate) %>% 
  mutate(day = str_extract(biosample1, "\\d{3}"))

recipients_day29_families <- recipient_ASVs %>%
  filter(day == "029", replicate == curr_replicate) %>% #  add replicate 1 filter
  select(biosample1, Family, OTU, relAbundance) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  distinct()

recipients_day36_families <- recipient_ASVs %>%
  filter(day == "036", replicate == curr_replicate) %>% #  add replicate 1 filter
  select(biosample1, Family, OTU, relAbundance) %>%
  mutate(subject1 = str_sub(biosample1, 1, -5)) %>% 
  distinct()

family_abundance_29 <- recipients_day29_families %>% 
group_by(subject1, Family) %>% 
  summarize(fam_abundance_29 = sum(relAbundance))


family_abundance_36 <- recipients_day36_families %>% 
  group_by(subject1, Family) %>% 
  summarize(fam_abundance_36 = sum(relAbundance))


lost_family_abundance <- family_abundance_29 %>% 
  left_join(family_abundance_36, by = c("subject1", "Family")) %>% 
  mutate(fully_lost = ifelse(is.na(fam_abundance_36), 1, 0))

fully_lost_families <- lost_family_abundance %>% 
  filter(fully_lost == 1)


fully_lost_family_plot <- lost_family_abundance %>% 
  select(subject1, Family, fam_abundance_29, fully_lost) %>% 
  filter(fully_lost == 1) %>% 
  ggplot(aes(x = fct_reorder(Family, - fam_abundance_29), y = fam_abundance_29, fill = subject1)) +
  geom_col(color = "black")+
  labs(x = "Family",
       y = "Day 29 relAbundance of lost family")+
  facet_wrap(~subject1)+
  scale_fill_manual(values = PALETTE.SUBJECT)+
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 5, angle = 90),
        legend.position = "none")


savePNGPDF(paste0(OUTDIR, "fully_lost_fams_day36"), fully_lost_family_plot , 3, 3)


