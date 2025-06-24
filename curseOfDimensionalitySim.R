library(randomForest)
library(tidyverse)
library(MASS)
library(furrr)
library(rdist)
library(progressr)

# First set up the simulation
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

sampleSizes <- c(50,100,150)
features <- c(10,50,100,150,200,250)
blockSize <- 10
repeats <- 1:10
minkowskiP <- c(0.25,0.5,1,2,1000)
rho <- 0.75
perms <- 25

params <- expand_grid(sampleSizes,features,blockSize,repeats,minkowskiP,rho,perms) %>% 
  mutate(seed = row_number())

plan(multisession, workers = 10)
future::nbrOfWorkers()


runSim <- function(samples,features,seed,blockSize,rho,run,perms,minkowski,p){
  
  set.seed(seed)
  
  # Generate Data + Prep
  
  data <- blocks(samples,features,blockSize,rho)
  
  
  data <- cbind(data,rowSums(data)) # outcome
  data <- cbind(data,runif(samples)) # baseline feature
  
  colnames(data) <- c(paste0("x",1:features),"y","B")
  data <- data.frame(data)
  
  learn_idx <- sample(1:nrow(data), round(0.75*nrow(data)))
  
  train <- data[learn_idx, ]
  test <- data[-learn_idx, ]
  
  # Fit Forest
  forest <- randomForest(y ~ .,data = train, ntree = samples, keep.forest = T, keep.inbag = T)
  
  # Compute Distances for one correlated feature and baseline feature
  toWritePVals <- matrix(0,ncol = 1, nrow = perms)
  
  for (permIdx in 1:perms){
    toPerm <- test[,!names(test) %in% c("y")]
    avgDists <- matrix(0,ncol = 2, nrow = nrow(toPerm))
    
    # Compute average leaf community distance upon permuting for correlated feature (1st WLOG)
    toPerm[,1] <- sample(toPerm[,1], replace = F)
    
    leafCommunity <- rep(0,nrow(toPerm))
    totalDist <- rep(0,nrow(toPerm))
    
    
    for (i in 1:samples){ # For each tree...
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
        trainDataInLeaf <- train[trainingIndices[bootstrapIndicesInLeaf],!names(train) %in% c("y")]
        
        testDataInLeaf <- toPerm[which(leafIndicesTest == leafIdx),]
        
        if (minkowski == 1000){
          pairwiseDistances <- cdist(trainDataInLeaf, testDataInLeaf, metric = "chebyshev")
        } else {
          pairwiseDistances <- cdist(trainDataInLeaf, testDataInLeaf, metric = "minkowski", p = minkowski)
        }
        
        
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
    
    avgDists[,1] <- totalDist/leafCommunity
    
    toPerm <- test[,!names(test) %in% c("y")]
    
    # Compute average leaf community distance upon permuting for correlated feature (1st WLOG)
    toPerm[,features+1] <- sample(toPerm[,features+1], replace = F)
    
    leafCommunity <- rep(0,nrow(toPerm))
    totalDist <- rep(0,nrow(toPerm))
    
    for (i in 1:samples){ # For each tree...
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
        trainDataInLeaf <- train[trainingIndices[bootstrapIndicesInLeaf],!names(train) %in% c("y")]
        
        testDataInLeaf <- toPerm[which(leafIndicesTest == leafIdx),]
        
        if (minkowski == 1000){
          pairwiseDistances <- cdist(trainDataInLeaf, testDataInLeaf, metric = "chebyshev")
        } else {
          pairwiseDistances <- cdist(trainDataInLeaf, testDataInLeaf, metric = "minkowski", p = minkowski)
        }
        
        
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
    
    avgDists[,2] <- totalDist/leafCommunity
    
    # Now compute p-values
    resamples <- 10000
    
    diffs <- avgDists[,1] - avgDists[,2]
    testStat <- mean(diffs)
    distribution <- rep(31,resamples)
    for (r in 1:resamples){
      signs <- sample(c(-1,1),size = nrow(avgDists), replace = T)
      simDiffs <- signs*diffs
      distribution[r] <- mean(simDiffs)
    }
    
    moreExtreme <- sum(abs(distribution) >= abs(testStat))
    pVal <- (moreExtreme+1)/(resamples + 1)
    toWritePVals[permIdx,1] <- pVal
    
  }

  toWrite <- data.frame(toWritePVals)
  colnames(toWrite) <- "pVal"
  toWrite %>% 
    mutate(minkowski = ifelse(minkowski == 1000,'chebyshev',paste('minkowski',minkowski)),
           samples = samples, features = features,
           blockSize = blockSize,run = run) -> toWrite
  # Write to File
  filenamePVals <- paste0("/Users/aaronfoote/COURSES/Spurs2014/Thesis Research/bigSim/CoD2/pVals-",samples,"-",features,"-",blockSize,"-",minkowski,"-",run,'.csv')
  write.csv(toWrite, filenamePVals, row.names = FALSE)
  
  p()
}

with_progress({
  p <- progressor(steps = nrow(params))
  result <- future_pmap(
    list(
      params[['sampleSizes']],
      params[['features']],
      params[['seed']],
      params[['blockSize']],
      params[['rho']],
      params[['repeats']],
      params[['perms']],
      params[['minkowskiP']]),
    runSim, p = p)
})

