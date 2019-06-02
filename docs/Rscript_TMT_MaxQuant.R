#################################################################
#################### read in the evidence data #################
#################################################################
library(MSstatsTMT)
library(tidyr)
library(dplyr)

# First, get protein ID information
proteinGroups <- read.table("data/data_MaxQuant_TMT/proteinGroups.txt", sep = "\t", header = TRUE)

# Read in MaxQuant file: evidence.txt
evi <- read.table("data/data_MaxQuant_TMT/evidence.txt", sep="\t", header=TRUE)
colnames(evi)

# unique Raw.file, which is TMT run.
unique(evi$Raw.file)

###### the expected annotation file
# Read in annotation including condition and biological replicates: MaxQuant_annotation.csv
annot.maxquant <- read.csv("data/data_MaxQuant_TMT/MaxQuant_annotation.csv", header = TRUE)
head(annot.maxquant)

###### prepare the annotation for each run
runs <- unique(evi$Raw.file) # MS runs
Run_info <- data.frame(Run = runs) # initialize the run file 

## Add the mixture and technical replicate information
Run_info <- Run_info %>% 
  separate(Run, c("Mixture", "TechRepMixture"), sep="_0", remove = FALSE) 
Run_info

## clean the file and add fraction information 
Run_info <- Run_info %>% 
  mutate(Mixture = gsub("161117_SILAC_HeLa_UPS1_TMT10_", "", Mixture),
         TechRepMixture = gsub(".raw", "", TechRepMixture),
         Fraction = "F1")
Run_info

###### prepare the annotation for the channels in each mixture
colnames(evi)

channels <- c("channel.0", "channel.1", "channel.2", "channel.3", "channel.4", "channel.5", "channel.6", "channel.7", "channel.8", "channel.9")
mixtures <- unique(Run_info$Mixture)
mixtures

## create the file with channel information in each mixture
Group_info <- expand.grid(channels, mixtures)
colnames(Group_info) <- c("Channel", "Mixture")
head(Group_info)

## save the channel file and fill in the condition and biological replicate information manually
write.csv(Group_info, file = "Group_info.csv", row.names = FALSE)

## Now the condition information should be available in the file
Group_info_filled <- read.csv(file = "data/data_MaxQuant_TMT/Group_info_mq.csv")
head(Group_info_filled)


###### combine run information and channel information 
annotation <- full_join(Run_info, Group_info_filled)
nrow(annotation)

head(annotation)

#################################################################
####### Preprocessing with MaxQtoMSstatsTMTFormat function ########
#################################################################

?MaxQtoMSstatsTMTFormat

# reformating and pre-processing for MaxQuant output.
input.maxquant <- MaxQtoMSstatsTMTFormat(evidence=evi, 
                                         annotation=annot.maxquant,
                                         proteinGroups=proteinGroups)
head(input.maxquant)

save(input.maxquant, file='data/data_MaxQuant_TMT/input.maxquant.rda')


#################################################################
########## Protein summarization and normalization ##############
#################################################################
load(file='data/data_MaxQuant_TMT/input.maxquant.rda')

?proteinSummarization

# channel level protein summarization and between-run normalization
quant.maxquant <- proteinSummarization(data = input.maxquant, 
                                 method = "msstats", 
                                 normalization = TRUE)

save(quant.maxquant, file='data/data_MaxQuant_TMT/quant.maxquant.rda')

head(quant.maxquant)

# Profile plot for the normalized data 
dataProcessPlotsTMT(data.psm=input.maxquant, # PSM-level data
                    data.summarization=quant.maxquant, # protein-level data
                    type='ProfilePlot', # choice of visualization
                    width = 15,
                    height = 5,
                    address="data/data_MaxQuant_TMT/mq_norm_") 


#################################################################
############ Statistical modeling and testing ###################
#################################################################
?groupComparisonTMT

# compare all the possible pairs of conditions
test.maxquant.pairwise <- groupComparisonTMT(data = quant.maxquant, 
                                       contrast.matrix = "pairwise",
                                       remove_norm_channel = TRUE, # remove norm channels
                                       moderated = TRUE, # do moderated t test
                                       adj.method = "BH") # multiple comparison adjustment

# show the comparisons
unique(test.maxquant.pairwise$Label)

# Show test result
# Label : which comparison is used
# log2FC : estimated log2 fold change between two conditions (the contrast)
# adj.pvalue : adjusted p value
head(test.maxquant.pairwise)

save(test.maxquant.pairwise, file='data/data_MaxQuant_TMT/test.maxquant.pairwise.rda')
write.csv(test.maxquant.pairwise, file='data/data_MaxQuant_TMT/test.maxquant.pairwise.csv')


head(test.maxquant.pairwise)

#### check the significant results ####
# select subset of rows with adj.pvalue < 0.05
SignificantProteins <- 
  test.maxquant.pairwise[test.maxquant.pairwise$adj.pvalue <= 0.05 ,]

SignificantProteins

nrow(SignificantProteins)

# select subset of rows with adj.pvalue < 0.05 and absoluate log2FC > 1
SignificantProteinsLargeFC <- SignificantProteins[abs(SignificantProteins$log2FC) > 1 ,]
SignificantProteinsLargeFC

nrow(SignificantProteinsLargeFC)
