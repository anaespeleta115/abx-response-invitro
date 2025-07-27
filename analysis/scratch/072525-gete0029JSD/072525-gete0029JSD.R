# Load data
source("C:/abx-response-invitro/analysis/plotDefaults.R")
source("C:/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")



beta_diversity <- read_delim("C:/abx-response-invitro/data/e0026_e0029_speciesBeta.txt.gz")
data(beta_diversity)


beta_diversity_wide <- beta_diversity %>%
  select(sample1, sample2, method, value) %>%
  pivot_wider(names_from = method, values_from = value)

beta_diversity_wide <- beta_diversity_wide %>%
  filter(!grepl("e0029", sample1), !grepl("e0029", sample2))



recipients <- c("XBA-029", "XBA-036", "XBA-064", "XDA-029", "XDA-036", "XDA-064", "XEA-029", "XEA-036", "XEA-064", "XKA-029", "XKA-036", "XKA-064")

e0026_filtered <- e0026 %>% 
  select(sample, biosample1) %>% 
  unique()

donor_communities <- e0029 %>%
  filter(biosample2 != "blank", 
         biosample2 != "B-mix",
         !grepl("\\+", biosample2)) %>%
  pull(biosample2) %>%
  unique()

# join sample1 metadata and check
e0026_e0029_beta_1 <- beta_diversity_wide %>%
  left_join(e0026_filtered, by = c("sample1" = "sample")) %>% 
  filter(biosample1 %in% recipients)

# if the columns exist, rename them
e0026_e0029_beta_1 <- e0026_e0029_beta_1 %>%
  dplyr::rename(
    component1 = biosample1
  )

# Join sample2 metadata
e0026_e0029_beta_final <- e0026_e0029_beta_1 %>%
  left_join(e0026_filtered, by = c("sample2" = "sample")) %>% 
  filter(biosample1 %in% donor_communities)


# Rename columns to have biosample1 and biosample
e0026_e0029_beta_final <- e0026_e0029_beta_final %>%
  dplyr::rename(
    biosample1 = component1,
    biosample2 = biosample1
  )




extract_passage <- function(x) {
  ifelse(
    str_detect(x, "-mix-"),
    NA_character_,
    str_extract(x, "-[A-Z]-[08]-[A-Z0-9]+") %>%
      str_extract("[08]")
  )
}

# Add passage extraction
e0026_e0029_beta_final <- e0026_e0029_beta_final %>%
  mutate(
    passage1 = extract_passage(sample1),
    passage2 = extract_passage(sample2)
  ) %>% 
  filter(passage1 == 8, passage1 == passage2)

e0026_e0029_beta_final <- e0026_e0029_beta_final %>% 
  select(-passage1, -passage2, -bray, -jaccard) %>% 
  mutate(mixture = paste(biosample2, biosample1, sep = "+"))







