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
library(ggcorrplot)


load("/Users/aaronfoote/Downloads/Singh.rda")

xTrain <- Singh$X
yTrain <- Singh$y
xTest <- Singh$Xt
yTest <- Singh$yt

# Select Features with top 5% variance
variances <- data.frame(apply(xTrain,2,var))
names(variances) <- "colVariance"
variances %>% mutate(colIdx = row_number()) %>% arrange(desc(colVariance)) %>%
  head(0.05*nrow(variances)) %>% pull(colIdx) -> highVarColIdxs

smallerTrain <- xTrain[,highVarColIdxs]
smallerTest <- xTest[,highVarColIdxs]

# Apply PCA to identify number of components we want
pca <- prcomp(smallerTrain)

var_explained_df <- data.frame(PC=1:min(nrow(smallerTrain),ncol(smallerTrain)),
                               var_explained=(pca$sdev)^2/sum((pca$sdev)^2))

var_explained_df %>% head(15) %>% 
  ggplot(aes(x=PC,y=var_explained, group = 1)) + geom_point(size = 4) + geom_line(linewidth = 2) + 
  labs(y = "Proportion of Variance Explained", x = "Component") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2))

ggsave(filename = '../../../RVisuals/SinghScree.tiff',device='tiff', dpi=1200)


# We want 7! Now we want to identify best values of λ_1 and λ_2 that provide a good balance between sparsity and explained variance

xTrain <- Singh$X

# Select Features with top 5% variance
variances <- data.frame(apply(xTrain,2,var))
names(variances) <- "colVariance"
variances %>% mutate(colIdx = row_number()) %>% arrange(desc(colVariance)) %>%
  head(0.05*nrow(variances)) %>% pull(colIdx) -> highVarColIdxs

smallerTrain <- xTrain[,highVarColIdxs]

l1s <- c(10^(seq(-7, 2)),5*10^(seq(-7, 2)))
l2s <- c(1e-6,1e-5,1e-4,1e-3,1e-2)

params <- expand_grid(l1s,l2s)

trySPCA <- function(l1,l2){
  spcaOutput <- elasticnet::spca(x = smallerTrain, K = 7, sparse = "penalty", type = "predictor",
                                 para = rep(l1,7), lambda = l2, max.iter = 2500)
  sum(spcaOutput$loadings > 0)
  sum(spcaOutput$pev)
  
  lineToWrite <- matrix(c(l1,l2,sum(spcaOutput$pev),sum(spcaOutput$loadings > 0)), nrow = 1)
  colnames(lineToWrite) <- c("l1","l2","propVarExplained","numNonZeroWeights")
  write.csv(
    data.frame(lineToWrite), 
    paste0("/Users/aaronfoote/COURSES/Spurs2014/Thesis Research/bigSim/spcaSim/Singh/",l1,"-",l2,".csv"),
    row.names = F)
}


future_pmap(list(params[['l1s']],params[['l2s']]),trySPCA)




# From examination of the visual produced in visuals.R, we find the values below to perform best.

spcaOutput <- spca(smallerTrain, K = 7, type = 'predictor', sparse = 'penalty', 
                   lambda = 1e-4, para = rep(3e-5,7), max.iter = 5000)

print(sum(spcaOutput$pev))



trainReduced <- smallerTrain %*% spcaOutput$loadings
testReduced <- smallerTest %*% spcaOutput$loadings

trainReduced <- cbind(trainReduced,factor(yTrain))
trainReduced <- cbind(trainReduced,runif(nrow(trainReduced)))
testReduced <- cbind(testReduced,factor(yTest))
testReduced <- cbind(testReduced,runif(nrow(testReduced)))

colnames(trainReduced) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","y","B")
colnames(testReduced) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","y","B")
trainReduced <- data.frame(trainReduced)
testReduced <- data.frame(testReduced)

forest <- randomForest(factor(y) ~ .,data = trainReduced, keep.forest = T,
                       keep.inbag = T, importance = T)
conditionalImportances <- permimp(forest, xdata = testReduced,
                                  y = 'y', conditional = T, nperm = 25, do_check = F)

importances <- tibble(conditionalImportances$values,
                      forest$importance[,3],
                      c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","B"))
colnames(importances) <- c("Conditional","Unconditional","Feature")

importances %>% 
  pivot_longer(names_to = "Method",cols = c(Conditional,Unconditional)) %>%
  ggplot(aes(x = value, y = Feature, fill = Method)) + geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(values = c("Conditional" = "darkgreen","Unconditional" = "orangered2")) +
  labs(x = "Average decrease in accuracy") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2),
        text = element_text(family = "serif"))
ggsave(filename = './COURSES/Spurs2014/IJCAI Paper/Tiff Visuals/SinghImps.tiff',device='tiff', dpi=1200)


perms <- 25
toWritePVals <- matrix(0,ncol = ncol(testReduced) - 1, nrow = perms)

for (permIdx in 1:perms){
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
    diffs <- avgDists[,i] - avgDists[,8]
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


colnames(toWritePVals) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","B")
data.frame(toWritePVals) %>% 
  pivot_longer(cols = PC1:B,names_to = "Feature", values_to = "p-value") %>%
  ggplot(aes(x = `p-value`, y = Feature, fill = Feature)) + geom_boxplot() + 
  xlim(c(0,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2),
        text = element_text(family = "serif"))
ggsave(filename = './COURSES/Spurs2014/IJCAI Paper/Tiff Visuals/SinghPVals.tiff',device='tiff', dpi=1200)

ggcorrplot(round(cor(testReduced %>% dplyr::select(-y)),2),
           lab = T,outline.color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA,linewidth = 2),
        text = element_text(family = "serif"))
ggsave(filename = './COURSES/Spurs2014/IJCAI Paper/Tiff Visuals/SinghCorr.tiff',device='tiff', dpi=1200)