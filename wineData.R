
# Load in packages
library(tidyverse)
library(ggcorrplot)
library(randomForest)
library(permimp)
library(rdist)

# Read in the Data
wine <- read_csv("/Users/aaronfoote/Downloads/wine/wine.data", col_names = FALSE) 
colnames(wine) <- c("Class","Alcohol","Malicacid","Ash","AshAlcalinity","Magnesium",
                    "TotalPhenols","Flavanoids","NonflavanoidPhenols","Proanthocyanins",
                    "ColorIntensity","Hue","ProteinConcentration","Proline")
wine$Class <- as.factor(wine$Class)

# Just for some preliminary exploration, let's take a look at the pairwise correlations between the variables.
wineCorrelation <- cor(wine %>% select(-Class))

ggcorrplot(round(wineCorrelation,2), lab = T,outline.color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2))

# Train-Test Split
set.seed(31)
wine['B'] <- runif(nrow(wine))
sample <- sample.int(n = nrow(wine), size = floor(.75*nrow(wine)), replace = F)
train <- wine[sample, ]
test  <- wine[-sample, ]


perms <- 25
treeCount <- 200
forest <- randomForest(Class ~ .,data = train, keep.forest = T,
                       keep.inbag = T, importance = T, ntree = treeCount)

# Calculate the conditional PFI scores (Strobl et al., 2008 Conditional variable importance for random forests)
conditionalImportances <- permimp(forest, xdata = test,
                                  y = 'Class', conditional = T, nperm = perms, do_check = F)

importances <- tibble(conditionalImportances$values,
                      forest$importance[,1],
                      c("Alcohol","Malicacid","Ash","Ash Alcalinity","Magnesium",
                        "Total Phenols","Flavonoids","Nonflavanoid Phenols","Proanthocyanins",
                        "Color Intensity","Hue","Protein Concentration","Proline","B"))
colnames(importances) <- c("Conditional","Unconditional","Feature")

# Combine PFI and cPFI scores for comparison
importances %>% 
  pivot_longer(names_to = "Method",cols = c("Conditional","Unconditional"),values_to = "Importance") %>%
  ggplot(aes(x = Importance, y = Feature)) + geom_bar(aes(fill = Method),stat = 'identity', position="dodge") +  
  labs(x = "Average decrease in accuracy") + 
  scale_fill_manual(values = c("Conditional" = "darkgreen","Unconditional" = "orangered2")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2),
        text = element_text(family = "serif"))

# I find that .tiff files retain very high resolution when converted to .png
ggsave(filename = '/Users/aaronfoote/COURSES/Spurs2014/IJCAI Paper/Tiff Visuals/wineImps.tiff',device='tiff', dpi=1200)


testSmaller <- test %>% select(-Class)
toWritePVals <- matrix(0,ncol = ncol(test) - 1, nrow = perms)

for (permIdx in 1:perms){
  toPerm <- test[,!names(test) %in% c("Class")]
  avgDists <- matrix(0,ncol = ncol(toPerm), nrow = nrow(toPerm))
  for (colIdx in 1:ncol(toPerm)){ # For each column...
    # Shuffle the column
    toPerm[,colIdx] <- sample_n(toPerm[,colIdx],size = nrow(toPerm), replace = F)
    
    leafCommunity <- rep(0,nrow(toPerm))
    totalDist <- rep(0,nrow(toPerm))
    
    for (i in 1:treeCount){ # For each tree...
      # Identify the training observations used to fit
      trainingIndices <- unname(which(forest$inbag[,i] > 0))
      # Get the frequency of each one
      trainingFrequencies <- unname(forest$inbag[trainingIndices,i])
      # Get the actual training data
      trainingObs <- train[trainingIndices,]
      # Identify which leaves this training subset ends up in
      leafIndicesBoot <- unname(attr(predict(forest, newdata = trainingObs, nodes = T), "nodes")[,i])
      
      # Identify which leaves the permuted data end up in
      leafIndicesTest <- unname(attr(predict(forest, newdata = toPerm, nodes = T), "nodes")[,i])
      
      uniqueLeaves <- unique(leafIndicesTest)
      
      for (leafIdx in uniqueLeaves){
        bootstrapIndicesInLeaf <- which(leafIndicesBoot == leafIdx)
        trainDataInLeaf <- train[trainingIndices[bootstrapIndicesInLeaf],!names(train) %in% c("Class")]
        
        testDataInLeaf <- toPerm[which(leafIndicesTest == leafIdx), !names(toPerm) %in% c("Class")]
        
        pairwiseDistances <- cdist(trainDataInLeaf, testDataInLeaf, metric = "euclidean", p = 2)
        
        for (testDist in 1:ncol(pairwiseDistances)){
          totalDist[which(leafIndicesTest == leafIdx)[testDist]] <-
            totalDist[which(leafIndicesTest == leafIdx)[testDist]] +
            sum(pairwiseDistances[,testDist]*trainingFrequencies[bootstrapIndicesInLeaf])
          leafCommunity[which(leafIndicesTest == leafIdx)[testDist]] <-
            leafCommunity[which(leafIndicesTest == leafIdx)[testDist]] +
            sum(trainingFrequencies[bootstrapIndicesInLeaf])
        }
      }
    }
    
    avgDists[,colIdx] <- totalDist/leafCommunity
    
    # Don't forget to change it back!
    toPerm[,colIdx] <- testSmaller[,colIdx]
  }
  
  resamples <- 10000
  n <- nrow(avgDists)
  for (i in 1:ncol(avgDists)){
    diffs <- avgDists[,i] - avgDists[,ncol(test)-1]
    testStat <- mean(diffs)
    
    
    distribution <- rep(31,resamples)
    for (r in 1:resamples){
      signs <- sample(c(-1,1),size = nrow(avgDists), replace = T)
      simDiffs <- signs*diffs
      distribution[r] <- mean(simDiffs)
    }
    
    moreExtreme <- sum(abs(distribution) >= abs(testStat))
    pVal <- (moreExtreme+1)/(resamples + 1)
    toWritePVals[permIdx,i] <- pVal
  }
  
}

colnames(toWritePVals) <- c("Alcohol","Malicacid","Ash","Ash Alcalinity","Magnesium",
                            "Total Phenols","Flavonoids","Nonflavanoid Phenols","Proanthocyanins",
                            "Color Intensity","Hue","Protein Concentration","Proline","B")
data.frame(toWritePVals) %>% 
  pivot_longer(cols = Alcohol:B,names_to = "Feature", values_to = "p-value") %>%
  ggplot(aes(x = `p-value`, y = Feature, fill = Feature)) + geom_boxplot() + 
  xlim(c(0,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2),
        text = element_text(family = "serif"))

ggsave(filename = '/Users/aaronfoote/COURSES/Spurs2014/IJCAI Paper/Tiff Visuals/winePVals.tiff',device='tiff', dpi=1200)
