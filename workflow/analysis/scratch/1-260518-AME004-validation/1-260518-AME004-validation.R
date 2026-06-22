source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-loadeACE010eAME004/1-260406-loadeACE010eAME004.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260428-geteAME004_ASV_lists/1-260428-geteAME004_ASV_lists.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260430-runPCA/1-260430-runPCA.R")

OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260518-AME004-validation/out/"


# Plot each taxa's dynamics for each subject:

trajectories <- recipient_ASVs %>% 
  ggplot(aes(x = day, y = relAbundance, fill = Family, group = Family, color = OTU)) +
  geom_area(position = "stack", color = "black") +
  facet_wrap(~subject) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Study day", y = "Relative abundance") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT

trajectories

# I think this plot only shows the OTUs that were present in all of the study timepoints


otu_trajectories <- recipient_ASVs %>%
  ggplot(aes(x = day, y = relAbundance, color = OTU, group = OTU)) +
  geom_line() +
  facet_wrap(~subject) +
  labs(x = "Study day", y = "Relative abundance") +
  ylim(0, 1)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT


savePNGPDF(paste0(OUTDIR, "otu_trajectories"), otu_trajectories, 5, 5)

# try on the ACE010 data

ACE_otu_trajectories <- eACE010_data %>% 
  filter(biosample == "XCA-028" | biosample %in% recipient_ASVs$biosample1) %>% 
  ggplot(aes(x = day, y = relAbundance, color = OTU, group = OTU)) +
  geom_line() +
  facet_wrap(~subject) +
  labs(x = "Study day", y = "Relative abundance") +
  ylim(0, 1)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT

ACE_otu_trajectories

savePNGPDF(paste0(OUTDIR, "ACE_otu_trajectories"), ACE_otu_trajectories, 5, 5)


# try correlating relative abundances

eACE010_filtered <- eACE010_data %>% filter(biosample == "XCA-028" | biosample %in% recipient_ASVs$biosample1)
recipient_filtered <- recipient_ASVs %>% select(-mixture, -biosample2, -replicate) %>% rename(biosample = biosample1)

colnames(eACE010_filtered)
colnames(recipient_filtered)

both_trajectories <- eACE010_filtered %>%  
  full_join(recipient_filtered, by = c("biosample", "Family", "OTU"))

correlate_trajectories <- both_trajectories %>% 
  ggplot(aes(x = relAbundance.x, y = relAbundance.y, color = Family)) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  labs(
    x = "eACE010 relative abundance",
    y = "eAME004 relative abundance"
  ) +
  scale_x_log10()+
  scale_y_log10()+
  coord_cartesian(xlim = c(1e-4, 1), ylim = c(0.001, 0.5)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  DEFAULTS.THEME_PRINT

correlate_trajectories
  
savePNGPDF(paste0(OUTDIR, "correlate_trajectories"), correlate_trajectories, 5, 5)

# Plot correlation between AME004 replicates





# Plot JSD for same mixes vs. different mixes in AME004


# Use JSD to compare same sample, different replicates, then compare different samples, different replicates.

## UNCOMMENT THIS ONLY WHEN NECESSARY TO RE-COMPUTE JSD

distMethods <- c("jsd")

calculateBeta <- function(data, distMethod) {
  # Calculate the distance matrix using the specified method.
  betaRaw <- distance(data, method=distMethod)


  beta <- as.matrix(betaRaw)
  beta <- as.data.frame(beta)
  beta$sample1 <- rownames(beta)

  beta <- beta %>%
    pivot_longer(-sample1, names_to="sample2", values_to="value")
  beta <- beta %>%
    filter(sample1 != sample2) %>%
    mutate(method=distMethod)
}

# Calculate the distance matrix for all of the specified methods on the species abundances.
# Combine the distance matrices for all methods.
betaSpecies <- do.call(rbind, lapply(distMethods, function(distMethod) {
  print(distMethod)
  calculateBeta(subset_obj, distMethod)
}))
# Export the distance matrix generated for all of the sample pairs
# using all of the specified methods.
write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta.txt.gz"))


# Add sample IDs back in

beta_diversity <- read_delim("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260518-AME004-validation/out/speciesBeta.txt.gz")

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

  







