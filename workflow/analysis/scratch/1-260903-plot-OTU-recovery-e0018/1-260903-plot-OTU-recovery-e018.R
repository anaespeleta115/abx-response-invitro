source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260903-plot-OTU-recovery-e0018/out/"


# read in csv file generated at the end of 1-260903-percent-recovered-e0018.R script:
recovered_otus <- read.csv("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260903-percent-recovered-e0018/out/all_households_recovered_OTUs.csv")


# use pivot_longer to create a new column (gained_lost_OTUs) and another column

recovered_otus_longer <- recovered_otus |> 
  select(household, lost_otus, recovered_otus, recovery_percentage) |> 
  mutate(household = fct_reorder(household, recovery_percentage, .desc = FALSE)) |> 
  pivot_longer(
    cols = c(lost_otus, recovered_otus), # The columns you want to collapse
    names_to = "Gain/Loss",  # Name for the new column of 
    values_to = "Num.OTUs"   # Name for the column containing the data values
  )

# plot two separate points for each household in one plot, one for the lost OTUs and one for the recovered OTUs
p_recovered_otus_longer <- recovered_otus_longer %>% 
  ggplot(aes(x = household, y = Num.OTUs, color = `Gain/Loss`))+
  geom_point()+
  geom_line(aes(group=`Gain/Loss`))+
  labs(
    x = "Household",
    y = "Number of OTUs"
  )+
  ylim(0, 50)+
  DEFAULTS.THEME_PRINT

p_recovered_otus_longer

savePNGPDF(paste0(OUTDIR, "p_recovered_otus_longer"), p_recovered_otus_longer, 3, 4)


# Try plotting the number of OTUs on the y axis and the two discrete variables lost and recovered on the x axis. Then connect them with a line. 
# Each line would be a different household

p_recovered_otus_bipartite <- recovered_otus_longer |> 
  ggplot(aes(x = as.factor(`Gain/Loss`), y = Num.OTUs, color = household)) +
  geom_point()+
  geom_line(aes(group=household))+
  labs(
    x = "",
    y = "Number of OTUs"
  ) +
  scale_x_discrete(
    labels = c(
      "lost_otus" = "Lost",
      "recovered_otus" = "Recovered"
    )
  ) +
  ylim(0, 50)+
DEFAULTS.THEME_PRINT

p_recovered_otus_bipartite

savePNGPDF(paste0(OUTDIR, "p_recovered_otus_bipartite"), p_recovered_otus_bipartite, 3, 3.5)


