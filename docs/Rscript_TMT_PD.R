#################################################################
#################### read in the psm level data #################
#################################################################
library(MSstatsTMT)
library(tidyr)
library(dplyr)

raw.pd <- read.delim("data/data_ProteomeDiscoverer_TMT/spikedin_PSMs.txt")
colnames(raw.pd) # check the colnames

# total number of unique protein name
proteins <- unique(raw.pd$Protein.Accessions)
length(proteins)

# show the spiked-in proteins
proteins[grepl("ups",proteins)]

# total number of unique peptide names
length(unique(raw.pd$Annotated.Sequence))

# unique Spectrum.File, which is TMT run.
unique(raw.pd$Spectrum.File)


###### the expected annotation file
annot.pd <- read.csv(file="data/data_ProteomeDiscoverer_TMT/PD_Annotation.csv")

###### prepare the annotation for each run
runs <- unique(raw.pd$Spectrum.File) # MS runs
Run_info <- data.frame(Run = runs) # initialize the run file 
?separate

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
colnames(raw.pd)

channels <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")
mixtures <- unique(Run_info$Mixture)
mixtures

## create the file with channel information in each mixture
Group_info <- expand.grid(channels, mixtures)
colnames(Group_info) <- c("Channel", "Mixture")
head(Group_info)

## save the channel file and fill in the condition and biological replicate information manually
write.csv(Group_info, file = "Group_info.csv", row.names = FALSE)

## Now the condition information should be available in the file
Group_info_filled <- read.csv(file = "data/data_ProteomeDiscoverer_TMT/Group_info_pd.csv")
head(Group_info_filled)


###### combine run information and channel information 
annotation <- full_join(Run_info, Group_info_filled)
nrow(annotation)

head(annotation)

#################################################################
####### Preprocessing with PDtoMSstatsTMTFormat function ########
#################################################################

?PDtoMSstatsTMTFormat

# reformating and pre-processing for PD output.
processed.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation=annotation)
head(processed.pd)

# total number of proteins
proteins <- unique(processed.pd$ProteinName)
length(proteins)

# show the spiked-in proteins
proteins[grepl("ups",proteins)]

save(processed.pd, file='data/data_ProteomeDiscoverer_TMT/processed.pd.rda')
write.csv(processed.pd, file='data/data_ProteomeDiscoverer_TMT/processed.pd.csv', row.names = FALSE)


#################################################################
########## Protein summarization and normalization ##############
#################################################################
load(file='data/data_ProteomeDiscoverer_TMT/input.pd.rda')

?proteinSummarization

# channel level protein summarization and between-run normalization
quant.pd <- proteinSummarization(data = input.pd, 
                                 method = "msstats", 
                                 normalization = TRUE)

save(quant.pd, file='data/data_ProteomeDiscoverer_TMT/quant.pd.rda')

head(quant.pd)

# Profile plot is good visualization to check individual measurements.
# if you have many MS runs, adjust width of plot (make wider)

# Profile plot for the normalized data 
dataProcessPlotsTMT(data.psm=input.pd, # PSM-level data
                    data.summarization=quant.pd, # protein-level data
                    type='ProfilePlot', # choice of visualization
                    width = 15,
                    height = 5,
                    address="data/data_ProteomeDiscoverer_TMT/pd_norm_") 

# Instead of making all profile plots for all proteins, we can make plot for individual protein.
dataProcessPlotsTMT(data.psm=input.pd, # PSM-level data
                    data.summarization=quant.pd, # protein-level data
                    type='ProfilePlot', # choice of visualization
                    width = 15,
                    height = 5,
                    which.Protein = 'P35221',
                    address="data/data_ProteomeDiscoverer_TMT/pd_norm_P35221_") 

# channel level protein summarization without normalization
quant.pd.nonorm <-proteinSummarization(data = input.pd, 
                                       method = "msstats", 
                                       normalization = FALSE)

# Profile plot for the data without normalization
dataProcessPlotsTMT(data.psm = input.pd, # PSM-level data
                    data.summarization = quant.pd.nonorm, # protein-level data
                    type = 'ProfilePlot', # choice of visualization
                    width = 15,
                    height = 5,
                    originalPlot = FALSE,
                    which.Protein = 'P35221',
                    address="data/data_ProteomeDiscoverer_TMT/pd_noNorm_P35221_")

# Different summarization option (median summarization)
quant.pd.median <-proteinSummarization(data = input.pd, 
                                       method = "Median", 
                                       normalization = TRUE)

dataProcessPlotsTMT(data.psm=input.pd, # PSM-level data
                    data.summarization=quant.pd.median, # protein-level data
                    type='ProfilePlot', # choice of visualization
                    width = 15,
                    height = 5,
                    originalPlot = FALSE,
                    which.Protein = 'P35221',
                    address="data/data_ProteomeDiscoverer_TMT/pd_median_P35221_") 


#################################################################
############ Statistical modeling and testing ###################
#################################################################
?groupComparisonTMT

# compare all the possible pairs of conditions
test.pd.pairwise <- groupComparisonTMT(data = quant.pd, 
                                       contrast.matrix = "pairwise",
                                       remove_norm_channel = TRUE, # remove norm channels
                                       moderated = TRUE, # do moderated t test
                                       adj.method = "BH") # multiple comparison adjustment

# show the comparisons
unique(test.pd.pairwise$Label)

# Show test result
# Label : which comparison is used
# log2FC : estimated log2 fold change between two conditions (the contrast)
# adj.pvalue : adjusted p value
head(test.pd.pairwise)

# If you have multiple groups, you can assign any group comparisons you are interested in.
# check unique conditions and check order of condition information
# In this case, four different concentrations
unique(quant.pd$Condition)

# 'Norm' will be removed during tesing and should be not considered in the contrast
comparison1<-matrix(c(-1,0,0,1),nrow=1) # 0.5-0.125
comparison2<-matrix(c(0,-1,1,0),nrow=1) # 0.667-0.5
comparison<-rbind(comparison1, comparison2)
# Set the column names
colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
# Set the names of each row
row.names(comparison)<-c("1-0.125","0.667-0.5")

comparison

test.pd <- groupComparisonTMT(data = quant.pd, 
                              contrast.matrix = comparison,
                              remove_norm_channel = TRUE, # remove norm channels
                              moderated = TRUE, # do moderated t test
                              adj.method = "BH") # multiple comparison adjustment

colnames(test.pd)

head(test.pd)

save(test.pd, file='data/data_ProteomeDiscoverer_TMT/pd.result.rda')
write.csv(test.pd, file='data/data_ProteomeDiscoverer_TMT/testResult_pd.csv')


head(test.pd)

#### check the significant results ####
# select subset of rows with adj.pvalue < 0.05
SignificantProteins <- 
  test.pd[test.pd$adj.pvalue <= 0.05 ,]

SignificantProteins

nrow(SignificantProteins)

# select subset of rows with adj.pvalue < 0.05 and absoluate log2FC > 1
SignificantProteinsLargeFC <- SignificantProteins[abs(SignificantProteins$log2FC) > 1 ,]
SignificantProteinsLargeFC

nrow(SignificantProteinsLargeFC)


