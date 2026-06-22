# Load data
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/scratch/102025-loade0040data/102025-loade0040data.R")
source("~/Documents/GitHub/abx-response-invitro/workflow/analysis/plotDefaults.R")

e0041 <- fread("~/Documents/GitHub/bottom-up-community-mixtures/data/ps_all_e0041.txt.gz")

e0041_obj <- readRDS("~/Documents/GitHub/bottom-up-community-mixtures/data/ps_all_e0041.rds")


# Define color palette
my_colors <- readRDS("~/Documents/GitHub/abx-response-invitro/data/familyColorPalette.rds") 

# Define replicate used for colonzation function
curr_replicate <- 1

# set limit of detection
limit_of_detection <- 0.001

# assign community label and replicate based on wells
e0041 <- e0041 %>% 
  mutate(
    experiment = case_when(
    plate == "e0041-A-5"~ "XEA_exp",
    plate == "e0041-B-5"~ "XBA_exp",
    TRUE ~ NA_character_
  ),
    recipient = case_when(
    str_detect(well, "A") ~ "blank",
    str_detect(well, "B|C") ~ "pre-abx",
    str_detect(well, "D|E") ~ "post-abx-V1",
    str_detect(well, "F|G") ~ "post-abx-V2",
    str_detect(well, "H") ~ "blank",
    TRUE ~ NA_character_
  ),
  replicate = case_when(
    str_detect(well, "A|B|D|F") ~ "1",
    str_detect(well, "C|E|G|H") ~ "2"
  ),
  donor = case_when(
    str_detect(well, "10$") ~ "XKB-029",
    str_detect(well, "11$") ~ "super",
    str_detect(well, "12$") ~ "blank",
    str_detect(well, "1$") ~ "XBB-029",
    str_detect(well, "2$") ~ "XCB-029",
    str_detect(well, "3") ~ "XDB-029",
    str_detect(well, "4") ~ "XEB-029",
    str_detect(well, "5") ~ "XFB-029",
    str_detect(well, "6") ~ "XGB-029",
    str_detect(well, "7") ~ "XHB-029",
    str_detect(well, "8") ~ "XIB-029",
    str_detect(well, "9") ~ "XJB-029",
    TRUE ~ NA_character_
  ),
  mixture = paste(recipient, donor, sep = "+")) 
  
# filter for subject XEA and XBA
e0041.A <- e0041 %>% filter(experiment == "XEA_exp", relAbundance>limit_of_detection) 
e0041.B <- e0041 %>% filter(experiment == "XBA_exp") 

# single donor, single recipients
e0041_control_donors <- e0041.B %>% filter(recipient == "blank", donor != "blank", replicate == 1, relAbundance>limit_of_detection) # replicate 2 donors look better
e0041_control_recipients <- e0041.A %>% filter(donor == "blank", recipient != "blank", replicate == 1, relAbundance>limit_of_detection) # replicate 1 recipients looks better

# calculate richness for each group of communities
e0041_control_donors_richness <- e0041_control_donors %>% group_by(replicate, mixture, well) %>% summarize(richness = sum(relAbundance>0.001))
e0041_control_recipients_richness <- e0041_control_recipients %>% group_by(replicate, mixture, well) %>% summarize(richness = sum(relAbundance>0.001))
e0041_mixture_richness_XEA <- e0041.A %>% group_by(replicate, mixture, well) %>% summarize(richness=sum(relAbundance>0.001))


# get all mixtures
mixture_ids <- unique(e0041.A %>% 
                      filter(recipient != "blank", donor != "blank") %>% 
                      pull(mixture))

# get ordered wells
e0041.A <- e0041.A %>%
  mutate(well = sprintf("%s%02d",
                        substr(well, 1, 1),
                        as.integer(substring(well, 2))))

rows <- LETTERS[1:8]
cols <- sprintf("%02d", 1:12)
wells <- as.vector(sapply(rows, function(r) paste0(r, cols)))
