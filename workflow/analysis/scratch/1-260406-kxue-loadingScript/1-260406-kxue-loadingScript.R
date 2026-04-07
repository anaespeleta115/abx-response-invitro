library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library(patchwork)
library("phyloseq")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/background/background.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")


OUTDIR <- "~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/1-260406-kxue-loadingScript/out/"

# Import phyloseq object for new communities.
data <- read.table("~/Documents/GitHub/abx-response-invitro/data/eACE010-eAME004ps_all.txt.gz",
                   header=TRUE, stringsAsFactors = FALSE)

# Extract community data.
data_eACE010 <- data %>%
  filter(round2plate=="A", round1index %in% c(0,1))
# Label based on row and column.
data_eACE010 <- data_eACE010 %>%
  mutate(well=substr(filename, 7, nchar(filename))) %>%
  mutate(well=substr(well, 1, nchar(well)-3))

# Calculate the diversity in each well.
data_eACE010_diversity <- data_eACE010 %>%
  group_by(metadata) %>%
  filter(relAbundance>0.001) %>%
  summarize(nASVs=n())

# Plot the diversity in each community.
data_eACE010_diversity %>%
  ggplot() +
  stat_ecdf(aes(x=nASVs)) +
  xlab("Number of ASVs") + ylab("Proportion of communities")

# Plot a stacked bar plot of composition in a select set of communities.
data_eACE010 <- data_eACE010 %>%
  filter(metadata %in% c("XBA-022", "XBA-029", "XBA-036", "XBA-064"),
         relAbundance>0.001) %>%
  ggplot() +
  geom_bar(aes(x=metadata, y=relAbundance, fill=Family), stat="identity", color="black")

savePNGPDF(paste0(OUTDIR, "kxue_script_eACE10_compositions"), data_eACE010, 8, 8)
