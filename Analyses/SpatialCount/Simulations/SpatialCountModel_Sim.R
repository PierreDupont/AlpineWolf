rm(list=ls())
library(nimble)


## ------ CUSTOM FUNCTIONS -----
dpoisSC_normal <- nimbleFunction(
  run = function( x = double(1),
                  lambda0 = double(0),
                  sigma = double(0),
                  s = double(2),
                  trapCoords = double(2),
                  oper = double(1, default = 1.0),
                  indicator = double(1),
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(0))
    
    ## Identify dimensions
    numIndividuals <- nimDim(s)[1]
    numTraps <- nimDim(trapCoords)[1]
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(numTraps)
    
    ## Loop over individuals
    for(i in 1:numIndividuals){
      if(indicator[i]>0){
        ## Loop over traps 
        for(j in 1:numTraps){
          d2 <- pow(trapCoords[j,1] - s[i,1], 2) + pow(trapCoords[j,2] - s[i,2], 2)
          lambda[j] <- lambda[j] + lambda0 * exp(alpha * d2) 
        }#if
      }#j
    }#i
    
    ## Calculate the log-probability of the vector of detections
    logProb <- 0.0 
    ## Loop over traps 
    for(j in 1:numTraps){
      logProb <- logProb + dpois( x = x[j],
                                  lambda = lambda[j]*oper[j],
                                  log = TRUE)
    }#j
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })

rpoisSC_normal <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  lambda0 = double(0),
                  sigma = double(0),
                  s = double(2),
                  trapCoords = double(2),
                  oper = double(1, default = 1.0),
                  indicator = double(1)
  ){
    ## Specify return type
    returnType(double(1))
    if(n!=1){print("rpoisSC_normal only allows n = 1; using n = 1")}
    
    ## Identify dimensions
    numIndividuals <- nimDim(s)[1]
    numTraps <- nimDim(trapCoords)[1]
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric( length = numTraps)
    
    ## Loop over individuals
    for(i in 1:numIndividuals){
      if(indicator[i]>0){
        ## Loop over traps 
        for(j in 1:numTraps){
          d2 <- pow(trapCoords[j,1] - s[i,1], 2) + pow(trapCoords[j,2] - s[i,2], 2)
          lambda[j] <- lambda[j] + lambda0 * exp(alpha * d2) 
        }#if
      }#j
    }#i
    
    ## Calculate the log-probability of the vector of detections
    out <- nimNumeric(numTraps)
    ## Loop over traps 
    for(j in 1:numTraps){
      out[j] <- rpois(n = 1,
                      lambda = lambda[j]*oper[j])
      
    }#j
    
    ## Output
    return(out)
  })


## -----------------------------------------------------------------------------
## ----- 1. DATA SIMULATION ------
## -----   1.1. SIMULATION FUNCTION ------

# adapted sim function from oSCR, from Cat Sun
# changed to also output activity centers
sim.SC <- function( n.individuals = 100,                ## Number of groups
                    n.occasions = 1,                ## Number of detection occasions
                    lambda0 = 0.75,             ## Baseline probability of visit 
                    sigma = 0.5,           ## Scale parameter
                    traps_dim = c(10,10),    ## Detector grid size
                    trapCoords = NULL, 
                    buffer = NULL){
  ##-- Set up trap locations
  if(is.null(trapCoords)) {
    trapCoords <- expand.grid(X = seq(1, traps_dim[1], by = 1), 
                              Y = seq(1, traps_dim[2], by = 1))
  }
  n.traps <- nrow(trapCoords)
  
  ##-- Calculate buffer size
  if (is.null(buffer)) { buffer <- sigma * 3}
  
  ##-- Set up habitat boundaries of the rectangular habitat
  Xl <- min(trapCoords[, 1] - buffer)
  Xu <- max(trapCoords[, 1] + buffer)
  Yl <- min(trapCoords[, 2] - buffer)
  Yu <- max(trapCoords[, 2] + buffer)
  
  ##-- Sample group activity centers
  sx <- runif(n.individuals, Xl, Xu)
  sy <- runif(n.individuals, Yl, Yu)
  S <- cbind(sx, sy)
  
  ##-- Calculate average density in the population
  Dens <- n.individuals/((Xu - Xl) * (Yu - Yl))
  
  ##-- Function to calculate distances
  e2dist <- function (x, y){
    if (!is.matrix(x)) 
      x <- as.matrix(x)
    if (!is.matrix(y)) 
      y <- as.matrix(y)
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  ##-- Calculate distance between all traps and activity centers
  D <- e2dist(S, trapCoords)
  
  ##-- Calculate alpha for faster calculation
  alpha1 <- 1/(2 * sigma * sigma)
  
  ##-- Calculate group- and trap-specific visit probabilities
  lambda <- lambda0 * exp(-alpha1 * D * D)
  
  ##-- Set-up observation matrix y:
  ##-- the dimensions are n.traps*K and it contains the number of wolves 
  ##-- photographed for each trap and occasion
  Y <- array(NA, c(n.individuals,n.traps,n.occasions))
  dimnames(Y) <- list("groups" = c(1:n.individuals),
                      "traps" = paste("trap", 1:n.traps),
                      "occasions" = 1:n.occasions)
  ##-- Loop over occasions
  for(i in 1:n.individuals){
    for(j in 1:n.traps){
      for(k in 1:n.occasions){
        ##-- Sample number of wolves detections 
        Y[i,j,k] <- rpois(n = 1,lambda = lambda[i,j])
      }#k
    }#j
  }#i
  ##-- Output
  return(list( "Y" = apply(Y,c(2,3),sum),
               "n.individuals" = n.individuals,
               "n.traps" = n.traps,
               "n.occasions" = n.occasions,
               "S" = S,
               "trapCoords" = trapCoords,
               "habitatBoundaries" = list("x" = c(Xl,Xu),
                                          "y" = c(Yl,Yu))))
}



## -----  1.2. SIMULATE DATA ------
testSim <- sim.SC( n.individuals = 100,                ## Number of groups
                   lambda0 = 0.75,             ## Baseline probability of visit 
                   sigma = 0.5,           ## Scale parameter
                   traps_dim = c(5,10))



## -----  1.3. PREPARE NIMBLE DATA ------
##-- Nimble data
nimData <- list( y = c(testSim$Y))

##-- Nimble constants
nimConstants <- list( M = testSim$n.individuals*2,
                      J = testSim$n.traps,
                      oper = rep(1,testSim$n.traps),
                      trapCoords = as.matrix(testSim$trapCoords),
                      xlim = testSim$habitatBoundaries$x,
                      ylim =  testSim$habitatBoundaries$y)

##-- Nimble inits
z.init <- c( rep(1,testSim$n.individuals),
             rep(0,testSim$n.individuals))

nimInits <- list( s = rbind(testSim$S,
                            testSim$S),
                  lambda0 = runif(1,0,3),
                  psi = 0.5,
                  sigma = 0.5,          
                  z = z.init)

##-- Nimble parameters
params <- c("N", "lambda0", "psi", "sigma")



## -----------------------------------------------------------------------------
## ----- 2. MODEL FITTING ------
## -----  2.1. SPATIAL COUNT MODEL V.1 -----
## This is the basic formulation which loops over all individuals and detectors
SC_model <- nimbleCode({
  ##---- SPATIAL PROCESS 
  for(i in 1:M){ 
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
  }#g
  
  ##---- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)
  for(i in 1:M){ 
    z[i] ~ dbern(psi)
  }#i 								
  
  ##---- DETECTION PROCESS 
  sigma ~ dunif(0,10)
  alpha <- -1/(2*sigma*sigma)
  lambda0 ~ dunif(0,10)
  
  for(j in 1:J) {
    for(i in 1:M){ 
      distsq[i,j] <- (s[i,1] - trapCoords[j,1])^2 + (s[i,2] - trapCoords[j,2])^2
      lambda[i,j] <- lambda0 * exp(alpha*distsq[i,j]) * z[i] 
    }#i
    bigLambda[j] <- sum(lambda[1:M,j]) 
    y[j] ~ dpois(bigLambda[j]*oper[j])
  }#j
  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
})


model <- nimbleModel( code = SC_model,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
system.time(Cmodel$calculate())
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
system.time(out1 <- runMCMC(CmodelMCMC,
                            niter = 500,
                            nchains = 2,
                            samplesAsCodaMCMC = T))




## -----  2.2. SPATIAL COUNT MODEL V.2 -----
## This version uses a nimble function to speed-up calculations on an individual 
## basis. If the individual is not availabel for detection, calculation of 
## distances to detectors and detection rates are skipped.
calculateLambda <- nimbleFunction(
  run = function( lambda0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  indicator = double(0), 
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(1))
    
    ## Identify dimensions
    numTraps <- nimDim(trapCoords)[1]
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(length = numTraps)
    
    ## Shortcut if individual is not available for detection
    if (indicator < 1){return(lambda)}
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma) 
    
    ## Loop over traps
    for(j in 1:numTraps){
      d2 <- pow(trapCoords[j,1] - s[1], 2) + pow(trapCoords[j,2] - s[2], 2)
      lambda[j] <- lambda0 * exp(alpha * d2)
    }
    return(lambda)
  })


SC_model_2 <- nimbleCode({
  ##---- SPATIAL PROCESS 
  for(i in 1:M){ 
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
  }#g
  
  ##---- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)
  for(i in 1:M){ 
    z[i] ~ dbern(psi)
  }#i 								
  
  ##---- DETECTION PROCESS 
  sigma ~ dunif(0,10)
  lambda0 ~ dunif(0,10)
  
  for(i in 1:M){ 
    lambda[i,1:J] <- calculateLambda(
      lambda0 = lambda0,
      sigma = sigma,
      s = s[i,1:2],
      trapCoords = trapCoords[1:J,1:2],
      indicator = z[i])
  }#i
  
  for(j in 1:J){
    bigLambda[j] <- sum(lambda[1:M,j]) 
    y[j] ~ dpois(bigLambda[j]*oper[j])
  }#j

  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
})

model2 <- nimbleModel( code = SC_model_2,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits)
Cmodel2 <- compileNimble(model2)
system.time(Cmodel2$calculate())
modelConf2 <- configureMCMC(model2)
modelConf2$addMonitors(params)
modelMCMC2 <- buildMCMC(modelConf2)
CmodelMCMC2 <- compileNimble(modelMCMC2, project = model2)
system.time(out2 <- runMCMC(CmodelMCMC2,
                            niter = 500,
                            nchains = 2,
                            samplesAsCodaMCMC = T))




## -----  2.3. SPATIAL COUNT MODEL V.3 -----
## This version uses a nimble distribution to calculate the log-prob of all the 
## detections at all the traps at the same time. 
dpoisSC_normal <- nimbleFunction(
  run = function( x = double(1),
                  lambda0 = double(0),
                  sigma = double(0),
                  s = double(2),
                  trapCoords = double(2),
                  oper = double(1, default = 1.0),
                  indicator = double(1),
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(0))
    
    ## Identify dimensions
    numIndividuals <- nimDim(s)[1]
    numTraps <- nimDim(trapCoords)[1]
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(numTraps)
    
    ## Loop over individuals
    for(i in 1:numIndividuals){
      if(indicator[i]>0){
        ## Loop over traps 
        for(j in 1:numTraps){
          d2 <- pow(trapCoords[j,1] - s[i,1], 2) + pow(trapCoords[j,2] - s[i,2], 2)
          lambda[j] <- lambda[j] + lambda0 * exp(alpha * d2) 
        }#if
      }#j
    }#i
    
    ## Calculate the log-probability of the vector of detections
    logProb <- 0.0 
    ## Loop over traps 
    for(j in 1:numTraps){
      logProb <- logProb + dpois( x = x[j],
                                  lambda = lambda[j]*oper[j],
                                  log = TRUE)
    }#j
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


SC_model_3 <- nimbleCode({
  ##---- SPATIAL PROCESS 
  for(i in 1:M){ 
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
  }#g
  
  ##---- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)
  for(i in 1:M){ 
    z[i] ~ dbern(psi)
  }#i 								
  
  ##---- DETECTION PROCESS 
  sigma ~ dunif(0,10)
  lambda0 ~ dunif(0,10)
  
  y[1:J] ~ dpoisSC_normal(
    lambda0 = lambda0,
    sigma = sigma,
    s = s[1:M,1:2],
    trapCoords = trapCoords[1:J,1:2],
    oper = oper[1:J],
    indicator = z[1:M])
  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
})


model3 <- nimbleModel( code = SC_model_3,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits)
Cmodel3 <- compileNimble(model3)
system.time(Cmodel3$calculate())
modelConf3 <- configureMCMC(model3)
modelConf3$addMonitors(params)
modelMCMC3 <- buildMCMC(modelConf3)
CmodelMCMC3 <- compileNimble(modelMCMC3, project = model3)
system.time(out3 <- runMCMC(CmodelMCMC3,
                            niter = 500,
                            nchains = 2,
                            samplesAsCodaMCMC = T))




## -----------------------------------------------------------------------------
