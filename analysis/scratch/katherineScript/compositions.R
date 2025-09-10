ibrary(tidyverse)
library(cowplot)
library(data.table)
library(forcats)
library(foreach)

theme_set(theme_cowplot())

# Output directory for plots.
outdir <- "analysis/220207-communityComposition/"
# Import plot default aesthetics.
source("config/plotDefaults.R")

# Import read summary statistics.
dataRaw <- read.table("DADA2_output/L5_summary.txt", header=TRUE,
                      stringsAsFactors = FALSE)
taxa <- colnames(dataRaw)
dataRaw <- rownames_to_column(dataRaw) %>% dplyr::rename(sample=rowname)

# Import sample metadata.
meta <- read.table("config/220118-e0023-e0026-e0029-e0030-16S-samplesheet.txt",
                   header=TRUE, stringsAsFactors = FALSE) %>%
  filter(group=="e0026-e0029-e0030")

# Import taxa color palette.
palette <- read.table("config/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE)

# Parse sample names.
dataRaw <- dataRaw %>%
  separate(sample, into=c("filename","stem"), sep="_")
dataRaw <- left_join(dataRaw, meta, by=c("filename"))
dataRaw <- dataRaw %>% mutate(row=substr(well,1,1), col=substr(well,2,length(well)))
dataRaw$col <- as.numeric(dataRaw$col)

# Tidy data.
data <- dataRaw %>% gather(key="taxa",value="count",taxa)

# Calculate relative abundance.
data <- data %>%
  group_by(sample) %>%
  mutate(relAbundance=count/sum(count),
         relAbundance=ifelse(is.na(relAbundance),0,relAbundance))

# Create short version of taxa names.
data <- data %>%
  mutate(taxaShort=gsub(".*\\.","",taxa))

# Bind colors.
data <- left_join(data, 
                  palette %>% dplyr::select(taxa, hex) %>% filter(taxa %in% data$taxa), 
                  by=c("taxa"))
# Extract color palette. Convert unknown taxa to dark gray.
taxaPalette <- data %>% ungroup() %>%
  dplyr::select(taxa, taxaShort, hex) %>%
  unique() %>% mutate(hex=ifelse(is.na(hex),"#615c5c", hex)) %>%
  arrange(taxa)
taxaPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaPaletteList) <- taxaPalette$taxa
# Create a version of the color palette with shortened taxa names.
taxaPalette <- taxaPalette %>%
  mutate(taxaShort=gsub(".*\\.","",taxa))
taxaShortPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaShortPaletteList) <- taxaPalette$taxaShort

# Extract list of more prevalent taxa.
abundantTaxa <- (data %>% ungroup() %>% group_by(taxa) %>%
                   summarize(maxRelAbundance=max(relAbundance)) %>%
                   filter(maxRelAbundance>0.03))$taxa

# Plot key of abundant taxa.
p_paletteShort <- taxaPalette %>%
  filter(taxa %in% abundantTaxa) %>%
  ggplot() +
  geom_tile(aes(x=taxaShort, y=0, fill=factor(taxaShort))) +
  scale_fill_manual(values=taxaShortPaletteList) +
  guides(fill=FALSE) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(outdir,"palette-taxaShort.png"),
          p_paletteShort)  


# Plot community composition for e0026 communities ------------------------

# Plot stacked bar plots showing the relative abundance of families
# in the final communities.
p_passage8stacks <- data %>%
  filter(experiment=="e0026", passage==8) %>%
  separate(biosample1, into=c("subject","timepoint"), sep="-") %>%
  mutate(timepoint=as.numeric(timepoint)) %>%
  ggplot() +
  geom_bar(aes(x=factor(timepoint), y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~subject, scales="free", ncol=8) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Timepoint") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL
save_plot(paste0(outdir,"e0026-passage8-stacks.png"), p_passage8stacks,
          ncol=3, nrow=2)

# Plot stacked bar plots showing the relative abundance of families
# in the starting communities.
p_passage0stacks <- data %>%
  filter(experiment=="e0026", passage==0) %>%
  separate(biosample1, into=c("subject","timepoint"), sep="-") %>%
  mutate(timepoint=as.numeric(timepoint)) %>%
  ggplot() +
  geom_bar(aes(x=factor(timepoint), y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~subject, scales="free", ncol=8) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Timepoint") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL
save_plot(paste0(outdir,"e0026-passage0-stacks.png"), p_passage0stacks,
          ncol=3, nrow=2)

# Plot the change in family-level community composition over time
# for samples that were sequenced at every passage.
samplesEveryPassage <- data %>%
  filter(experiment=="e0026") %>%
  group_by(biosample1) %>%
  summarize(numPassages=n_distinct(passage)) %>%
  filter(numPassages==9)
p_communityStability <- data %>%
  filter(experiment=="e0026", 
         biosample1 %in% samplesEveryPassage$biosample1) %>%
  ggplot() +
  geom_bar(aes(x=passage, y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~biosample1, scales="free", ncol=4) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Passage") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL
save_plot(paste0(outdir,"e0026-stability-stacks.png"), p_communityStability,
          ncol=1.5, nrow=1)

p_communityStabilitySubset <- data %>%
  filter(experiment=="e0026", 
         biosample1 %in% samplesEveryPassage$biosample1,
         biosample1 %in% c("XGB-029","XJA-029")) %>%
  mutate(biosample1=ifelse(biosample1=="XGB-029","Sample 1","Sample 2")) %>%
  ggplot() +
  geom_bar(aes(x=passage, y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~biosample1, scales="free", ncol=1) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Passage") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL +
  theme(axis.title=element_text(size=6),
        axis.text=element_text(size=4.5),
        strip.text.x=element_text(size=6))
save_plot(paste0(outdir,"e0026-stability-stacks-large.png"), p_communityStabilitySubset,
          ncol=0.25, nrow=0.6)
abundantTaxaSelected <- (data %>%
                           filter(experiment=="e0026", 
                                  biosample1 %in% samplesEveryPassage$biosample1,
                                  biosample1 %in% c("XGB-029","XJA-029")) %>%
                           ungroup() %>% group_by(taxa) %>%
                           summarize(maxRelAbundance=max(relAbundance)) %>%
                           filter(maxRelAbundance>0.05))$taxa
# Plot key of abundant taxa in these two communities.
p_paletteShortSelected <- taxaPalette %>%
  filter(taxa %in% abundantTaxaSelected) %>%
  ggplot() +
  geom_tile(aes(x=taxaShort, y=0, fill=factor(taxaShort))) +
  scale_fill_manual(values=taxaShortPaletteList) +
  guides(fill=FALSE) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(outdir,"palette-taxaShort-selected.png"),
          p_paletteShortSelected, ncol=0.6, nrow=0.55) 

# Plot community composition for the e0029 communities --------------------

# Extract data for the e0029 mixture experiment.
datae0029 <- data %>%
  filter(experiment=="e0029", biosample1!="blank") %>%
  mutate(mixtureType=ifelse(biosample2=="blank","unmixed",
                            ifelse(grepl("[+]",biosample2),"double",
                                   ifelse(biosample2=="B-mix","all","mixed"))))
# Distinguish between the four replicates of the unmixed X[X]A-036 communities.
datae0029 <- datae0029 %>%
  mutate(replicate=ifelse(row %in% c("G","H"),replicate+2,replicate))

# Plot the composition of the donor communities.
# Extract the names of the donor communities.
donorCommunities <- unique((datae0029 %>%
                              filter(!grepl("[+]", biosample2), grepl("[X]", biosample2)))$biosample2)
# Display their compositions from the e0026 passaging experiment.
p_e0029stacksDonorCommunities <- data %>%
  filter((experiment=="e0026" & passage==8 &
            biosample1 %in% donorCommunities) |
           (experiment=="e0030" & passage==5 & biosample1=="B-mix")) %>%
  mutate(donorCommunity=paste(biosample1,replicate,sep="_")) %>%
  ggplot() +
  geom_bar(aes(x=factor(donorCommunity), y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Donor community") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))
save_plot(paste0(outdir,"e0029-donorCommunities-stacks.png"), 
          p_e0029stacksDonorCommunities,
          ncol=0.6, nrow=0.8)

# Plot the results of the initial community mixtures.
p_e0029stacksSingleMixtures <- datae0029 %>%
  mutate(biosample2=ifelse(biosample2=="B-mix","Bmix",biosample2)) %>%
  mutate(donorCommunity=paste(biosample2,replicate,sep="_")) %>%
  filter(mixtureType %in% c("unmixed","mixed","all")) %>%
  ggplot() +
  geom_bar(aes(x=factor(donorCommunity), y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~biosample1, scales="free", ncol=3) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Donor community") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))
save_plot(paste0(outdir,"e0029-singleMixtures-stacks.png"), 
          p_e0029stacksSingleMixtures,
          ncol=1.5, nrow=2)

# Reorder the naming of the double mixtures so that they fall
# between the two parent communities.
datae0029DoubleMixesRenamed <- datae0029 %>%
  mutate(biosample2=ifelse(biosample2=="B-mix","Bmix",biosample2)) %>%
  mutate(donorCommunity=paste(biosample2,replicate,sep="_")) %>%
  filter(biosample1 %in% c("XBA-036","XDA-036","XEA-036","XKA-036")) %>%
  separate(biosample2,into=c("biosample2A","biosample2B"),sep="[+]",remove=FALSE) %>%
  mutate(biosample2B=ifelse(is.na(biosample2B),"",biosample2B)) %>%
  mutate(biosample2=ifelse(grepl("[+]",biosample2),
                           paste(biosample2B, biosample2A,sep="+"),
                           biosample2)) %>%
  mutate(donorCommunity=paste(biosample2,replicate,sep="_"))

# Plot the results of the mixed community mixtures.
p_e0029stacksDoubleMixtures <- datae0029DoubleMixesRenamed %>%
  ggplot() +
  geom_bar(aes(x=factor(donorCommunity), y=relAbundance,
               fill=factor(taxa)), stat="identity") +
  facet_wrap(~biosample1, scales="free", ncol=2) +
  scale_fill_manual(values=taxaPaletteList) +
  guides(fill=FALSE) +
  xlab("Donor community") + ylab("Relative abundance") +
  DEFAULTS.THEME_ALL +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))
save_plot(paste0(outdir,"e0029-timepoint36Mixtures-stacks.png"), 
          p_e0029stacksDoubleMixtures,
          ncol=2, nrow=2)