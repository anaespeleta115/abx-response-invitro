source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-261113-get-eAME008.3-communities/1-261113-get-eAME008.3-communities.R")

# set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-261113-get-eAME008.3-compositions/out/"




# plot community compositions
p_bottom_up_communities <- bottomUp_data |> 
  ggplot(aes(x = community, y = `Relative Abundance`, fill = Family)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(title = "" , x = "Replicate 1 communities", y = "", fill = "Family") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT

p_bottom_up_communities

savePNGPDF(paste0(OUTDIR, "p_bottom_up_communities"), p_bottom_up_communities, 3, 3)


# add a new category for faceting: community type
bottomUp_data <- bottomUp_data |> 
  mutate(type = ifelse(community %in% c("A01", "A03", "A06"), "full-community", "knockout"))
         

# ------------------------ Loop through all communities to plot each one 

# create list of communities
communities <- unique(bottomUp_data |>  
                        pull(community))


# loop through all communities and plot each one
for (i in seq_along(communities)) {
  c <- communities[i]
  
  d <- bottomUp_data %>% 
    filter(.data$community == c)
  
  p <- ggplot(d, aes(x = community, y = `Relative Abundance`, fill = Family)) +
    geom_col(color = "black", linewidth = 0.1) +
    scale_fill_manual(values = my_colors) +
    theme_minimal() +
    labs(x = NULL)+
    theme(
      legend.position = "none"
    ) +
    DEFAULTS.THEME_PRINT
  
  
  savePNGPDF(paste0(OUTDIR, "p", c), p, 2, 2)
}


