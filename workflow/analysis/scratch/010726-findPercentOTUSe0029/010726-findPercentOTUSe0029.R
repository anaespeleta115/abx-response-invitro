# Load data
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072125-loade0029Data/072125-loade0029Data.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/plotDefaults.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getDifferentialColonizers/072325-getDifferentialColonizers.R")
source("/Users/aespelet/Documents/Github/abx-response-invitro/analysis/scratch/072325-getLostStrains/072325-getLostStrains.R")

# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/010726-findPercentOTUSe0029/out/"


# ------------------------------------------ %  differential colonizer OTUs out of the colonizer OTUs per family (day 36) -----------------------------------------------

colonizer_families <- actual_colonizers_results %>% 
  filter(colonized_day36 == 1, ) %>% 
  group_by(subject1, subject2, Family) %>% 
  summarize(fam_relAbundance = sum(relAbundance), 
            family_colonizer_otus = n()) %>% 
  filter(fam_relAbundance > 0.01)

diff_colonizer_families <- actual_colonizers_results %>% 
  filter(diff_colonizer_36 == 1) %>% 
  group_by(subject1, subject2, Family) %>% 
  summarize(family_diff_colonizer_otus = n())

percentage_differential <- colonizer_families %>% 
  left_join(diff_colonizer_families, by = c("subject1", "subject2", "Family")) %>% 
  mutate(family_diff_colonizer_otus = ifelse(is.na(family_diff_colonizer_otus), 0, family_diff_colonizer_otus), 
         percentage = as.numeric(family_diff_colonizer_otus) / as.numeric(family_colonizer_otus) * 100,
         above_50 = ifelse(percentage > 50, 1, 0))


# ------------------------------------------ plot -----------------------------------------------

# loop through all subject2's:

subject2_list <- c("XBB", "XCB", "XDB", "XEB", "XFB", "XGB", "XHB", "XIB", "XJB", "XKB")

for (s2 in subject2_list) {
  
  percent_colonizer_otus_family_plot <- percentage_differential %>% 
    filter(subject1 == "XKA", subject2 == s2) %>%
    ggplot(
      aes(
        x = fct_reorder(Family, -percentage),
        y = percentage,
        fill = Family,
        alpha = factor(above_50)
      )
    ) +
    geom_col(color = "black") +
    labs(
      x = "Family",
      y = "post-abx-only/post-abx colonizers XKA",
      fill = "Family"
    ) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 100)) +
    DEFAULTS.THEME_PRINT +
    facet_wrap(~subject2)+
    theme(
      axis.text.x = element_text(
        hjust = 1, vjust = 0.5, size = 6, angle = 90
      ),
      legend.position = "none"
    )
  
  savePNGPDF(
    paste0(OUTDIR, "percentage_diff_colonizer_otus_day36_XKA_", s2),
    percent_colonizer_otus_family_plot,
    3, 2
  )
}


