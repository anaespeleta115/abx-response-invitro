source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260812-loade0018data/1-260812-loade0018data.R")
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260818-plot-e0018-quality/out/"


# Here I plot the impurity percentage for each household using the two different samplesheets to corroborate that the later samplesheet
# has more inconsistencies than the later one. 

sept_samplesheet <- read.csv("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260817-calculate-e0018-quality/out/household_purity_20210901.csv")

feb_samplesheet <- read.csv("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260817-calculate-e0018-quality/out/household_purity_20220223.csv")



# ---------- Join both samplesheets by household

# rename columns to avoid redudancy
sept_samplesheet <- sept_samplesheet %>% 
  select(household, impurity_percentage) %>% 
  rename(impurity_percentage_sept = impurity_percentage)

feb_samplesheet <- feb_samplesheet %>% 
  select(household, impurity_percentage) %>% 
  rename(impurity_percentage_feb = impurity_percentage)


# join
compare_samplesheets <- sept_samplesheet %>% 
  left_join(feb_samplesheet, by = "household") %>% 
  mutate(a_to_j = ifelse(household %in% LETTERS[1:10], "A-J", "K-V"))



#  --------------- make a scatter plot, faceted by the two groups A-J and K-V
p_compare_samplesheets <- compare_samplesheets %>% 
  ggplot(aes(x = impurity_percentage_feb, y = impurity_percentage_sept))+
  geom_point()+
  labs(
    x = "Impurity Percentage February",
    y = "Impurity Percentage September"
  )+
  facet_wrap(~ a_to_j)

p_compare_samplesheets

savePNGPDF(paste0(OUTDIR, "p_compare_samplesheets"), p_compare_samplesheets, 3, 4.5)


# --------------- 

compare_samplesheets_longer <- compare_samplesheets %>% 
  pivot_longer(
  cols = c(impurity_percentage_feb, impurity_percentage_sept), # The columns you want to collapse
  names_to = "samplesheet_month",  # Name for the new column of 
  values_to = "impurity"   # Name for the column containing the data values
) 


p_compare_samplesheets_longer <- compare_samplesheets_longer %>% 
  ggplot(aes(x = samplesheet_month, y = impurity, fill = samplesheet_month))+
  geom_boxplot()+
  labs(
    x = "Samplesheet Month",
    y = "Impurity Percentage"
  )+
  facet_wrap(~ a_to_j)+
  DEFAULTS.THEME_PRINT+
  theme(
    legend.position = "none"
  )+
  scale_x_discrete(
    labels = c(
      "impurity_percentage_feb" = "Feb samplesheet",
      "impurity_percentage_sept" = "Sept samplesheet"
    )
  )

p_compare_samplesheets_longer

savePNGPDF(paste0(OUTDIR, "p_compare_samplesheets_longer"), p_compare_samplesheets_longer, 2, 3.5)





