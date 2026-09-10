source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260903-compare-JSD-recovery-e0018/out/"

# read in relevant comparisons dataset from JSD analysis on the autologous mixtures
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260825-get-e0018-beta-diversity/1-260825-get-e0018-beta-diversity.R")

# read in csv file generated at the end of 1-260903-percent-recovered-e0018.R script:
recovered_otus <- read.csv("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260903-percent-recovered-e0018/out/all_households_recovered_OTUs.csv")

# we have the recovery metric and the JSD metric. We want to try to correlate the two by household.


relevant_comparisons <- relevant_comparisons |> 
  filter(combination == "same_household_autologous", mix_type_2 == "A-pre to A-post-1") |> 
  mutate(household = str_sub(biosample1_1, 2, 2))
  
JSD_recovery <- relevant_comparisons |> 
  left_join(recovered_otus, by = "household") |> 
  mutate(household = fct_reorder(household, recovery_percentage, .desc = FALSE))

# plot two separate points for each household in one plot, one for the lost OTUs and one for the recovered OTUs
p_JSD_recovery <- JSD_recovery %>% 
  ggplot(aes(x = recovery_percentage, y = value, color = household))+
  geom_point()+
  geom_abline(stat = "identity")+
  labs(
    x = "Recovery %",
    y = "JSD"
  )+
  ylim(0, 0.5)+
  DEFAULTS.THEME_PRINT

p_JSD_recovery

savePNGPDF(paste0(OUTDIR, "p_JSD_recovery"), p_JSD_recovery, 3, 4)


