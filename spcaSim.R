library(elasticnet)
library(stats)
library(tidyverse)
library(MASS)
library(furrr)
library(progressr)
library(beepr)
library(randomForest)
library(permimp)
library(rdist)

blocks <- function(sampleSize,features,blockSize,rho){
  cov <- matrix(0,nrow = features,ncol = features)
  i <- 1
  while (i+blockSize-1 <= features) {
    cov[i:(i+blockSize-1),i:(i+blockSize-1)] <- rho
    i <- i + blockSize
  }
  
  for (i in 1:features){
    cov[i,i] <- 1
  }
  
  return(MASS::mvrnorm(n = sampleSize, mu = rep(0,features), Sigma = cov))
}
plan(multisession, workers = 10)
nbrOfWorkers()


spcaSimPlusTest <- function(trainSize,testSize,iteration,p){
  fullData <- blocks(trainSize + testSize,150,25,0.75)
  
  learn_idx <- sample(1:nrow(fullData), trainSize)
  
  train <- fullData[learn_idx, ]
  test <- fullData[-learn_idx, ]
  attempts <- 0
  l1 <- 5
  while (T){
    spcaOutput <- spca(train, K = 6, type = 'predictor', sparse = 'penalty', 
                       lambda = 5e-6, para = rep(l1,6), max.iter = 1000)
    goodColumns <- 0
    for (i in 1:6){
      if (sum(spcaOutput$loadings[,i] != 0) <= 25){goodColumns <- goodColumns + 1}
    }
    if (goodColumns == 6){
      break
    } else {
      if (attempts >= 12){
        p()
        return()
      }
      attempts <- attempts + 1
      l1 <- l1 + 0.5
    }
  }
  trainReduced <- train %*% spcaOutput$loadings
  testReduced <- test %*% spcaOutput$loadings
  
  # outcome
  yTrain <- rowSums(train) 
  yTest <- rowSums(test)
  
  trainReduced <- cbind(trainReduced,yTrain)
  trainReduced <- cbind(trainReduced,runif(trainSize))
  testReduced <- cbind(testReduced,yTest)
  testReduced <- cbind(testReduced,runif(testSize))
  
  colnames(trainReduced) <- c("PC1","PC2","PC3","PC4","PC5","PC6","y","B")
  colnames(testReduced) <- c("PC1","PC2","PC3","PC4","PC5","PC6","y","B")
  trainReduced <- data.frame(trainReduced)
  testReduced <- data.frame(testReduced)
  
  forest <- randomForest(y ~ .,data = trainReduced, keep.forest = T,
                         keep.inbag = T, importance = T, ntree = 100)
  conditionalImportances <- permimp(forest, xdata = testReduced,
                                    y = 'y', conditional = T, nperm = 25, do_check = F)
  
  importances <- tibble(conditionalImportances$values,
                        forest$importance[,1],
                        c("PC1","PC2","PC3","PC4","PC5","PC6","B"))
  colnames(importances) <- c("Conditional","Unconditional","Feature")
  
  filenameImps <- paste0("/Users/aaronfoote/COURSES/Spurs2014/Thesis Research/bigSim/spcaSim/apply/",iteration,"-imps",'.csv')
  write.csv(importances, filenameImps, row.names = FALSE)
  
  toWritePVals <- matrix(0,ncol = ncol(testReduced) - 1, nrow = 25)
  
  for (permIdx in 1:25){
    toPerm <- testReduced[,!names(testReduced) %in% c("y")]
    avgDists <- matrix(0,ncol = ncol(toPerm), nrow = nrow(toPerm))
    for (colIdx in 1:ncol(toPerm)){ # For each column...
      # Shuffle the column
      toPerm[,colIdx] <- sample(toPerm[,colIdx], replace = F)
      
      leafCommunity <- rep(0,nrow(toPerm))
      totalDist <- rep(0,nrow(toPerm))
      
      for (i in 1:100){ # For each tree...
        # Identify the training observations used to fit
        trainingIndices <- unname(which(forest$inbag[,i] > 0))
        # Get the frequency of each one
        trainingFrequencies <- unname(forest$inbag[trainingIndices,i])
        # Get the actual training data
        trainingObs <- trainReduced[trainingIndices,]
        # Identify which leaves this training subset ends up in
        leafIndicesBoot <- unname(attr(predict(forest, newdata = trainingObs, nodes = T), "nodes")[,i])
        
        # Identify which leaves the permuted data end up in
        leafIndicesTest <- unname(attr(predict(forest, newdata = toPerm, nodes = T), "nodes")[,i])
        
        uniqueLeaves <- unique(leafIndicesTest)
        
        for (leafIdx in uniqueLeaves){
          bootstrapIndicesInLeaf <- which(leafIndicesBoot == leafIdx)
          trainDataInLeaf <- trainReduced[trainingIndices[bootstrapIndicesInLeaf],!names(trainReduced) %in% c("y")]
          
          testDataInLeaf <- toPerm[which(leafIndicesTest == leafIdx), !names(toPerm) %in% c("y")]
          
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
      toPerm[,colIdx] <- testReduced[,colIdx]
    }
    
    resamples <- 10000
    n <- nrow(avgDists)
    for (i in 1:ncol(avgDists)){
      diffs <- avgDists[,i] - avgDists[,7]
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
  
  
  filenamePVals <- paste0("/Users/aaronfoote/COURSES/Spurs2014/Thesis Research/bigSim/spcaSim/apply/",iteration,"-pVals",'.csv')
  write.csv(toWritePVals, filenamePVals, row.names = FALSE)
  p()
}


trainSize <- 100
testSize <- 50
iterations <- 1:25
params <- expand_grid(trainSize,testSize,iterations)
with_progress({
  p <- progressor(steps = nrow(params))
  result <- future_pmap(
    list(
      params[['trainSize']],
      params[['testSize']],
      params[['iterations']]),
    spcaSimPlusTest, p = p)
})


# The code below simulates one instance of what is generated above, and is used to create the scree plots and bar charts in the thesis.

data <- blocks(150,150,25,0.75)

learn_idx <- sample(1:nrow(data), round(0.67*nrow(data)))

train <- data[learn_idx, ]
test <- data[-learn_idx, ]

pca <- prcomp(train)
spcaOutput <- spca(train, K = 6, type = 'predictor', sparse = 'penalty', 
                   lambda = 1e-4, para = rep(6,6), max.iter = 1000)

yTrain <- rowSums(train) # outcome
dimensionReducedData <- train %*% spcaOutput$loadings
dimensionReducedData <- cbind(dimensionReducedData,yTrain)
dimensionReducedData <- cbind(dimensionReducedData,runif(100))

colnames(dimensionReducedData) <- c("PC1","PC2","PC3","PC4","PC5","PC6","y","B")
dimensionReducedData <- data.frame(dimensionReducedData)


yTest <- rowSums(test)

transformedTestData <- test %*% spcaOutput$loadings
transformedTestData <- cbind(transformedTestData,yTest)
transformedTestData <- cbind(transformedTestData,runif(50))

colnames(transformedTestData) <- c("PC1","PC2","PC3","PC4","PC5","PC6","y","B")
transformedTestData <- data.frame(transformedTestData)

#### Visuals ----------------

var_explained_df <- data.frame(PC=as.character(1:100), var_explained=(pca$sdev)^2/sum((pca$sdev)^2))

var_explained_df %>% head(9) %>% 
  ggplot(aes(x=PC,y=var_explained, group = 1)) + geom_point(size = 4) + geom_line(linewidth = 2) + 
  labs(y = "Proportion of Variance Explained", x = "Component") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2))
ggsave(filename = '../../RVisuals/simScree.tiff',device='tiff', dpi=1200)

data.frame(pca$rotation[,1:6]) %>% 
  mutate(Feature = row_number()) %>% 
  pivot_longer(cols = PC1:PC6, names_to = "Component") %>%
  ggplot(aes(x = Feature, y = value)) + labs(y = "Weight") + 
  geom_bar(stat="identity") + facet_grid(rows = vars(Component)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 1),
        strip.background = element_rect(colour = 'black',fill = NA, linewidth = 1))
ggsave(filename = '../../RVisuals/pcaWeightsSim.tiff',device='tiff', dpi=1200)


data.frame(spcaOutput$loadings) %>% mutate(Feature = row_number()) %>%
  pivot_longer(cols = PC1:PC6, names_to = "Component") %>%
  ggplot(aes(x = Feature, y = value)) + labs(y = "Weight") + 
  geom_bar(stat="identity") + facet_grid(rows = vars(Component)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 1),
        strip.background = element_rect(colour = 'black',fill = NA, linewidth = 1))
ggsave(filename = '../../RVisuals/spcaWeightsSim.tiff',device='tiff', dpi=1200)