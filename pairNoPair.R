library(tidyverse)
library(MASS)
library(furrr)
library(rdist)
library(progressr)

blocks <- function(sampleSize,features,blockSize,rho){
  if (features < blockSize){ # This will happen when we have 2 features
    cov <- matrix(0,nrow = features,ncol = features)
    for (i in 1:features){
      cov[i,i] <- 1
    }
    cov[1,2] <- rho
    cov[2,1] <- rho
    
    return(MASS::mvrnorm(n = sampleSize, mu = rep(0,features), Sigma = cov))
    
  } else {
    cov <- matrix(0,nrow = features,ncol = features)
    i <- 1
    while (i+blockSize <= features) {
      cov[i:(i+blockSize-1),i:(i+blockSize-1)] <- rho
      i <- i + blockSize
    }
    
    for (i in 1:features){
      cov[i,i] <- 1
    }
    
    return(MASS::mvrnorm(n = sampleSize, mu = rep(0,features), Sigma = cov))
  }
}


plan(multisession, workers = 3)

sampleSizes <- c(10000)
features <- seq(2,10,2)
blockSize <- 3
rho <- c(0.25,0.5,0.75)
generators = c('standard MVN', 'uniform', 'blocks')

params <- expand_grid(sampleSizes,features,generators,rho,blockSize) %>% mutate(seed = row_number())
runRanking <- function(samples,features,seed,generator,blockSize,rho,p) {
  set.seed(seed)
  if (generator == 'standard MVN'){
    data <- mvrnorm(n = samples, mu = rep(0,features), Sigma = diag(ncol = features, nrow = features))
  }
  if (generator == 'uniform'){
    data <- as.matrix(mapply(runif,rep(samples,features)))
  }
  if (generator == 'blocks'){
    data <- blocks(samples,features,blockSize,rho)
  }
  
  permColumnIndices <- sample(ncol(data),2)
  permFirst <- data
  permFirst[,permColumnIndices[1]] <- sample(data[,permColumnIndices[1]])
  permSecond <- data
  permSecond[,permColumnIndices[2]] <- sample(data[,permColumnIndices[2]])
  
  distances <- cdist(permFirst,permSecond)
  
  freq_df <- as.data.frame(sapply(c(1:samples), (function (i) sum(distances[i,] <= distances[i,i]))))
  colnames(freq_df) <- c("Rank")
  filename <- paste0("/Users/aaronfoote/COURSES/Spurs2014/Thesis Research/simResults/",generator,"-",features,"-",samples,"-",rho,'.csv')
  write.csv(freq_df, filename, row.names = FALSE)
  
  p()
}



with_progress({
  p <- progressor(steps = nrow(params))
  result <- future_pmap(
    list(
      params[['sampleSizes']],params[['features']],
      params[['seed']],params[['generators']],
      params[['blockSize']],params[['rho']]),
    runRanking,
    p = p)
})