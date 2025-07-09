# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/041825-computeBetaPassages/out/"



### Define a function to extract passage information from raw sample ID

extract_passage <- function(x) {
  ifelse(
    str_detect(x, "-mix-"),
    NA_character_,
    str_extract(x, "-[A-Z]-[08]-[A-Z0-9]+") %>%
      str_extract("[08]")
  )
}


beta_diversity <- read_delim("C:/abx-response-invitro/analysis/scratch/041825-computeBetaDiversity/out/speciesBeta.txt.gz")

beta_diversity <- beta_diversity %>%
  select(sample1, sample2, method, value) %>%
  pivot_wider(names_from = method, values_from = value)

# ensure uniqueness in combined_data
combined_data_clean <- combined_day_data %>%
  distinct(sample, .keep_all = TRUE)

# join sample1 metadata and check
e0026_beta <- beta_diversity %>%
  left_join(combined_data_clean, by = c("sample1" = "sample"))

# if the columns exist, rename them
if (all(c("day", "subject", "antibiotic") %in% names(e0026_beta))) {
  e0026_beta <- e0026_beta %>%
    dplyr::rename(
      day1 = day,
      subject1 = subject,
      antibiotic1 = antibiotic
    )
} else {
  stop("Expected columns not found after join on sample1.")
}

# Join sample2 metadata
e0026_beta2 <- e0026_beta %>%
  left_join(combined_data_clean, by = c("sample2" = "sample"))


# Rename second set
if (all(c("day", "subject", "antibiotic") %in% names(e0026_beta2))) {
  e0026_beta2 <- e0026_beta2 %>%
    dplyr::rename(
      day2 = day,
      subject2 = subject,
      antibiotic2 = antibiotic
    )
} else {
  stop("Expected columns not found after join on sample2.")
}

# Add passage extraction
e0026_beta2 <- e0026_beta2 %>%
  mutate(
    passage1 = extract_passage(sample1),
    passage2 = extract_passage(sample2)
  )


### JSD between in vivo and in vitro communities (passage 0 vs. passage 8)


# Filter for values passage numbers and day 29 samples only

P0P8_JSD <- e0026_beta2 %>%
  filter(!is.na(passage1), !is.na(passage2), day1 == "029", day2 == "029")

# # Perform statistical test by looking at significance between passage pairs
# 
# left_join(P0P8_JSD %>%
#             mutate(combo=paste0(passage1,"_",passage2)) %>%
#             wilcox_test(jsd~combo),
#           P0P8_JSD %>%
#             mutate(combo=paste0(passage1,"_",passage2)) %>%
#             wilcox_effsize(jsd~combo),
#           by=c(".y.","group1","group2","n1","n2")) %>%
#   mutate(comparison=paste0(gsub("\n"," ",group1), " vs ", gsub("\n"," ",group2)),
#          summary=paste0("Wilcoxon rank-sum two-sided test, ", comparison, ": n=", n1+n2,
#                         ", r=", round(effsize,2), ", p=", p, ", adjusted p-value=", p.adj, ", ", p.adj.signif, ", custom adjusted p-value=", 4*p, " (number of groups: 4)"))

P0P8_JSD <- P0P8_JSD %>% 
  filter(passage1 == passage2)

p_P0P8_JSD <- ggplot(P0P8_JSD, aes(x = factor(passage1), y = jsd, fill = factor(passage1))) +
  geom_boxplot(width = 0.6, alpha = 0.85) +
  labs(
    title = "In Vivo vs In Vitro Community JSD (Within Passage)",
    x = "Passage Pair",
    y = "Jensen-Shannon Divergence",
    fill = "Passage"
  ) +
  theme(legend.position = "none")


savePNGPDF(paste0(OUTDIR, "betaByPassage"), p_P0P8_JSD, 6, 12)



