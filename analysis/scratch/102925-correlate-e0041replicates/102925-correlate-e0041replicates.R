# Load data
source("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102125-loade0041data/102125-loade0041data.R")


# Set output directory
OUTDIR <- "~/Documents/GitHub/abx-response-invitro/analysis/scratch/102925-correlate-e0041replicates/out/"


# Use JSD to compare same sample, different replicates, then compare different samples, different replicates.

### UNCOMMENT THIS ONLY WHEN NECESSARY TO RE-COMPUTE JSD
# 
# distMethods <- c("jsd")
# 
# calculateBeta <- function(data, distMethod) {
#   # Calculate the distance matrix using the specified method.
#   betaRaw <- distance(data, method=distMethod)
#   
# 
#   beta <- as.matrix(betaRaw)
#   beta <- as.data.frame(beta)
#   beta$sample1 <- rownames(beta)
# 
#   beta <- beta %>%
#     pivot_longer(-sample1, names_to="sample2", values_to="value")
#   beta <- beta %>%
#     filter(sample1 != sample2) %>%
#     mutate(method=distMethod)
# }
# 
# # Calculate the distance matrix for all of the specified methods on the species abundances.
# # Combine the distance matrices for all methods.
# betaSpecies <- do.call(rbind, lapply(distMethods, function(distMethod) {
#   print(distMethod)
#   calculateBeta(e0041_obj, distMethod)
# }))
# # Export the distance matrix generated for all of the sample pairs
# # using all of the specified methods.
# write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta.txt.gz"))


# Add sample IDs back in

beta_diversity <- read_delim("~/Documents/GitHub/abx-response-invitro/analysis/scratch/102925-correlate-e0041replicates/out/speciesBeta.txt.gz")

e0041_jsd <- beta_diversity %>%
  mutate(
    well1 = str_extract(sample1, "[A-H][0-9]{1,3}$"),
    well2 = str_extract(sample2, "[A-H][0-9]{1,3}$")
  )

# assign community labels and replicates based on wells
e0041_jsd <- e0041_jsd %>% 
  mutate(
    experiment1 = case_when(
      str_detect(sample1, "e0041-A-5") ~ "XEA_exp",
      str_detect(sample1, "e0041-B-5") ~ "XBA_exp",
      TRUE ~ NA_character_
    ),
    experiment2 = case_when(
      str_detect(sample2, "e0041-A-5") ~ "XEA_exp",
      str_detect(sample2, "e0041-B-5") ~ "XBA_exp",
      TRUE ~ NA_character_
    )) %>% 
  filter(experiment1 == "XEA_exp", experiment2 == "XEA_exp")


e0041_jsd_clean <- e0041_jsd %>% 
  mutate(
    recipient1 = case_when(
      str_detect(well1, "A") ~ "blank",
      str_detect(well1, "B|C") ~ "pre-abx",
      str_detect(well1, "D|E") ~ "post-abx-V1",
      str_detect(well1, "F|G") ~ "post-abx-V2",
      str_detect(well1, "H") ~ "blank",
      TRUE ~ NA_character_
    ),
    replicate_A = case_when(
      str_detect(well1, "A|B|D|F") ~ "1",
      str_detect(well1, "C|E|G|H") ~ "2"
    ),
    donor1 = case_when(
      str_detect(well1, "10$") ~ "XKB-029",
      str_detect(well1, "11$") ~ "super",
      str_detect(well1, "12$") ~ "blank",
      str_detect(well1, "1$") ~ "XBB-029",
      str_detect(well1, "2$") ~ "XCB-029",
      str_detect(well1, "3") ~ "XDB-029",
      str_detect(well1, "4") ~ "XEB-029",
      str_detect(well1, "5") ~ "XFB-029",
      str_detect(well1, "6") ~ "XGB-029",
      str_detect(well1, "7") ~ "XHB-029",
      str_detect(well1, "8") ~ "XIB-029",
      str_detect(well1, "9") ~ "XJB-029",
      TRUE ~ NA_character_
    ),
    mixture1 = paste(recipient1, donor1, sep = "+"
    ), 
    recipient2 = case_when(
      str_detect(well2, "A") ~ "blank",
      str_detect(well2, "B|C") ~ "pre-abx",
      str_detect(well2, "D|E") ~ "post-abx-V1",
      str_detect(well2, "F|G") ~ "post-abx-V2",
      str_detect(well2, "H") ~ "blank",
      TRUE ~ NA_character_
    ),
    replicate_B = case_when(
      str_detect(well2, "A|B|D|F") ~ "1",
      str_detect(well2, "C|E|G|H") ~ "2"
    ),
    donor2 = case_when(
      str_detect(well2, "10$") ~ "XKB-029",
      str_detect(well2, "11$") ~ "super",
      str_detect(well2, "12$") ~ "blank",
      str_detect(well2, "1$") ~ "XBB-029",
      str_detect(well2, "2$") ~ "XCB-029",
      str_detect(well2, "3") ~ "XDB-029",
      str_detect(well2, "4") ~ "XEB-029",
      str_detect(well2, "5") ~ "XFB-029",
      str_detect(well2, "6") ~ "XGB-029",
      str_detect(well2, "7") ~ "XHB-029",
      str_detect(well2, "8") ~ "XIB-029",
      str_detect(well2, "9") ~ "XJB-029",
      TRUE ~ NA_character_
    ),
    mixture2 = paste(recipient2, donor2, sep = "+")
    ) 


# Get plot distribution of same-subject, different replicates

same_mixes <- e0041_jsd_clean %>% 
  filter(mixture1 == mixture2) %>% 
  group_by(mixture1) %>%
  slice(1)

same_mixes_plot <- same_mixes %>% 
  ggplot(aes(x = 1, y = value)) +
  geom_boxplot(fill = "lightblue")  +
  labs(x = "", y = "Same-mix JSD") +
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "same_mixes_JSD_plot"), same_mixes_plot, 2.5, 2.5)

  

# Get distribution of different-subject, different replicates

diff_mixes <- e0041_jsd_clean %>%
  filter(mixture1 != mixture2) %>%
  group_by(mixture1) %>%
  slice(1)

diff_mixes_plot <- diff_mixes %>% 
  ggplot(aes(x = 1, y = value)) +
  geom_boxplot(fill = "lightgreen")  +
  labs(x = "", y = "Diff-mix JSD") +
  DEFAULTS.THEME_PRINT+
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())


savePNGPDF(paste0(OUTDIR, "diff_mixes_JSD_plot"), diff_mixes_plot, 2.5, 2.5)
