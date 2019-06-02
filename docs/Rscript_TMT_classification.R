#################################################################
#################### prepare the protein level data #############
#################################################################
library(tidyr)
library(dplyr)

# load data
load('data/data_ProteomeDiscoverer_TMT/quant.pd.rda')

# Pretend the two replicates within each condition and mixture are biological replicate
quant.pd <- quant.pd %>% 
  mutate(BioReplicate = paste(Mixture, Channel, sep="_"))


# protein quantification per subject
head(quant.pd)

# use technical replicate 2 and 3 as training data
quant.pd.per.subject <- quant.pd %>% filter(TechRepMixture != "1") %>% 
  group_by(Protein, BioReplicate) %>% 
  summarise(Abundance = median(Abundance, na.rm = TRUE)) %>%
  spread(BioReplicate, Abundance)

train_abun <- quant.pd.per.subject
colnames(train_abun)

# make protein abundance matrix
proteins <- train_abun$Protein
train_abun <- train_abun[, -1]
train_abun <- t(train_abun)
colnames(train_abun) <- proteins
dim(train_abun) # there are 50 rows (each row for subject) and 50 columns (one column per protein)

# get annotation information
colnames(quant.pd)
train_anno <- quant.pd %>% select(BioReplicate, Condition)
train_anno <- unique(train_anno)
train_anno <- as.data.frame(train_anno)
dim(train_anno) # there are 50 rows (each row for subject)
## remove the normalization channels
train_abun <- train_abun[train_anno$Condition != "Norm",]
train_anno <- train_anno[train_anno$Condition != "Norm",]

# Deal with missing values
# First option: remove the samples with missing values
dim(na.omit(train_abun))

# Second option: impute the missing values with miminal value
random.imp <- function (a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  # imputed[missing] <- 0 # with zero
  # imputed[missing] <- median(a.obs) # with median values
  imputed[missing] <- min(a.obs) # with minimal values
  return (imputed)
}

pMiss <- function(x){
  sum(is.na(x))/length(x)*100
}

# Only keep the subjects with less than 5% missing values
subjectmissing <- apply(train_abun, 1, pMiss)
train_abun <- train_abun[subjectmissing <= 5, ]
dim(train_abun)
# make sure the subject order in train_abun and train_anno consistent
train_anno <- train_anno[train_anno$BioReplicate %in% rownames(train_abun),] #remvoe the filtered subjects
train_abun <- train_abun[train_anno$BioReplicate,]

# Impute the missing values with minimum value per protein
imputed_train_abun <- apply(train_abun, 2, function(x) random.imp(x))
imputed_train_abun <- as.data.frame(imputed_train_abun)

sum(is.na(imputed_train_abun))


#################################################################
############# Principal components analysis (PCA) ###############
#################################################################
?prcomp
# rows are proteins and columns are subjects
pc <- prcomp(imputed_train_abun)

# Inspect PCA object
summary(pc)
names(pc)

# Check the proportion of explained variance
percent_var <- pc$sdev^2/sum(pc$sdev^2)
barplot(percent_var, xlab="Principle component", ylab="% of variance")

cum_var <- cumsum(pc$sdev^2/sum(pc$sdev^2))
barplot(cum_var, xlab="Principle component", ylab="Cumulative % of variance" )

# head(pc$x)
library(ggplot2)

# Let's visualize PC1 vs PC2 in scatterplot
ggplot(aes(x=PC1, y=PC2), data=data.frame(pc$x))+
  geom_point(size=4, alpha=0.5)+
  theme_bw()

# Create PC1 vs PC2 scatterplot with Condition colors
ggplot(aes(x=PC1, y=PC2, color=Condition), data=data.frame(pc$x, Condition=train_anno$Condition))+
  geom_point(size=4, alpha=0.5)+
  theme_bw()


#################################################################
################## Cluster analysis: Heatmap ####################
#################################################################
ht.data <- t(imputed_train_abun)
# check the class
class(ht.data)

# Change the font of row and column label
heatmap(ht.data, cexRow = 0.3, cexCol = 0.4)

library(marray)
my.colors <- c(maPalette(low = "darkblue", high = "white", k = 7)[-7],
               "white", 
               maPalette(low = "white", high = "darkred", k = 7)[-1])

heatmap(ht.data, cexRow = 0.3, cexCol = 0.4, col = my.colors)

# Don't do cluster on rows
heatmap(ht.data, cexRow = 0.3, cexCol = 0.4, col = my.colors, Rowv = NA)
# Don't do cluster on columns
heatmap(ht.data, cexRow = 0.3, cexCol = 0.4, col = my.colors, Colv = NA)

# Color bar for group information
unique(train_anno$Condition)
group.color <- rep("blue", nrow(imputed_train_abun))
group.color[train_anno$Condition == "1"] <- "red" 
group.color[train_anno$Condition == "0.667"] <- "yellow" 
group.color[train_anno$Condition == "0.5"] <- "orange" 
heatmap(ht.data, ColSideColors=group.color, cexRow = 0.3, cexCol = 0.4, Rowv = NA)


#################################################################
############# Classification and prediction #####################
#################################################################
# Set random seed to make results reproducible:
set.seed(430)
#install.packages("randomForest")
# Load library
library(randomForest)
?randomForest

# add group information to the training data
imputed_train_abun$Condition <- droplevels(train_anno$Condition)

# randomForest dosen't allow special symbol in the protein name
colnames(imputed_train_abun) <- gsub("-", "", colnames(imputed_train_abun))

# fit random forest
rf=randomForest(Condition ~ . , data = imputed_train_abun, importance=TRUE)
rf

# Prepare validation data 
# Use technical replicate 1 as validation data
valid_abun <- quant.pd %>% filter(TechRepMixture == "1") %>% 
  select(Protein, BioReplicate, Abundance) %>%
  spread(BioReplicate, Abundance)

proteins <- valid_abun$Protein
valid_abun <- valid_abun[, -1]
valid_abun <- t(valid_abun)
colnames(valid_abun) <- proteins
dim(valid_abun) 

# get annotation information
colnames(quant.pd)
valid_anno <- quant.pd %>% select(BioReplicate, Condition)
valid_anno <- unique(valid_anno)
valid_anno <- as.data.frame(valid_anno)
dim(valid_anno) # there are 50 rows (each row for subject)

valid_abun <- valid_abun[valid_anno$BioReplicate,]
## remove the normalization channels
valid_abun <- valid_abun[valid_anno$Condition != "Norm",]
valid_anno <- valid_anno[valid_anno$Condition != "Norm",]

imputed_valid_abun <- apply(valid_abun, 2, function(x) random.imp(x)) # impute missing values
imputed_valid_abun <- as.data.frame(imputed_valid_abun)
# randomForest dosen't allow special symbol in the protein name
colnames(imputed_valid_abun) <- gsub("-", "", colnames(imputed_valid_abun))

# prediction on validation set
valid_anno$Condition <- droplevels(valid_anno$Condition)
valid_pred <- predict(rf, imputed_valid_abun)

# Validation set assessment #1: looking at confusion matrix
table(data=valid_pred,
      reference=valid_anno$Condition)
# calculate the predictive accuracy
misClasificError <- mean(valid_pred != valid_anno$Condition)
print(paste('Accuracy',1-misClasificError))
