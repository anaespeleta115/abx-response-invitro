# Load data
source("C:/abx-response-invitro/analysis/scratch/040325-loadData/loadData.R")

library(coin)

source("C:/abx-response-invitro/analysis/plotDefaults.R")

# Set output directory
OUTDIR <- "C:/abx-response-invitro/analysis/scratch/041825-computeAlphaPrePostAbx/out/"

# Set palette
pal <- pnw_palette("Sunset2", 2, type = "discrete")

### Define a function to extract passage information from raw sample ID


extract_passage <- function(x) {
  ifelse(
    str_detect(x, "-mix-"),
    NA_character_,
    str_extract(x, "-[A-Z]-[08]-[A-Z0-9]+") %>%
      str_extract("[08]")
  )
}



### Plot change in Alpha Diversity over Time


# Read in alpha diversity file
alpha_diversity <- read_delim("C:/abx-response-invitro/analysis/scratch/041825-computeAlphaDiversity/out/alpha_diversity.txt")


# Only consider samples contained within combined data
valid_samples <- combined_day_data%>%
  distinct(sample) %>%
  pull(sample)

# Pivot to have each statistic as its own column
alpha_diversity <- alpha_diversity %>%
  filter(alphaStat %in% c("Shannon", "Chao1")) %>%
  select(sample, alphaStat, value) %>%
  pivot_wider(names_from = alphaStat, values_from = value)

# Add a column with altered sample ID
combined_data_alpha <- combined_day_data %>%
  distinct(sample, biosample1, day, subject, antibiotic)%>%
  mutate(xsample = paste0("X", sample))

# Left join be that sample ID
e0026_alpha <- combined_data_alpha %>%
  left_join(alpha_diversity, by = c("xsample" = "sample"))

# Add day, antibiotic, and passage information
e0026_alpha <- e0026_alpha %>%
  mutate(
    day = factor(day, levels = sort(unique(day))),
    antibiotic = factor(antibiotic, levels = c(0, 1), labels = c("No", "Yes")),
    passage = extract_passage(xsample)
  )

e0026_alpha <- e0026_alpha %>% 
  filter(passage == 8)


# Perform statistical test

valid_comparisons <- c("Yes_029 vs Yes_036", "Yes_029 vs Yes_064")


wilcoxon_results <- left_join(e0026_alpha %>%
            mutate(abxDay=paste0(antibiotic,"_",day)) %>%
            rstatix::wilcox_test(Shannon~abxDay),
          e0026_alpha %>%
            mutate(abxDay=paste0(antibiotic,"_",day)) %>%
            rstatix::wilcox_effsize(Shannon~abxDay),
          by=c(".y.","group1","group2","n1","n2")) %>%
  mutate(comparison=paste0(gsub("\n"," ",group1), " vs ", gsub("\n"," ",group2)),
         summary=paste0("Wilcoxon rank-sum two-sided test, ", comparison, ": n=", n1+n2,
                        ", r=", round(effsize,2), ", p=", p, ", adjusted p-value=", p.adj, ", ", p.adj.signif,
                        ", custom adjusted p-value=", 2*p, " (number of groups: 8")) %>%
  filter(comparison %in% valid_comparisons)



p_alpha_time <- ggplot(
  e0026_alpha,
  aes(x = day, y = Shannon, fill = antibiotic)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.9,
    size = 1
)  +
  # geom_jitter(
  #   color = "black",
  #   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
  #   shape = 21,
  #   stroke = 0.3,
  #   size = 1.5,
  #   alpha = 0.8
  # ) +
  labs(
    title = "Change in Alpha Diversity Pre and Post Abx",
    x = "Study Day",
    y = "Alpha Diversity",
    fill = "Antibiotic"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_fill_manual(values=pal)

savePNGPDF(paste0(OUTDIR, "alphaPrePostAbx"), p_alpha_time, 6, 12)
