
## The goal is to further develop the spatial count model (AKA the unmarked SCR) 
## based on camera trap data fro group-living species. The general idea is to 
## leverage the information contained in camera-trap pictures when multiple 
## individuals are detected together.

## The assumptions are:
## - all individuals detected together belong to the same group.
## - all individuals in a group share the same AC and home-range size.


## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
# library(nimbleSCR)





## -----------------------------------------------------------------------------
## ----- 1. DATA SIMULATION ------
# Adapted sim function from Catherine Sun
# changed to simulate group activity centers and group detections instead of indivdiuals
sim.SGM <- function( n.groups = 10,        ## Number of groups in the population
                    n.occasions = 20,      ## Number of occasions
                    p0 = 0.75,             ## Baseline probability of visit 
                    sigma = 0.5,           ## Scale parameter
                    lambda = 3,            ## Mean group size in the population
                    alpha = 0.75,          ## Cohesion of the group (probability to detect X ids together)
                    traps_dim = c(5,5),    ## Detector grid size
                    trapCoords = NULL, 
                    buffer = NULL){
  require(extraDistr)
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
  sx <- runif(n.groups, Xl, Xu)
  sy <- runif(n.groups, Yl, Yu)
  S <- cbind(sx, sy)
  
  ##-- Sample group sizes
  GS <- rtpois( n = n.groups,
                lambda = lambda,
                a = 0)
  #plot(S,cex = GS)
  
  ##-- Calculate average group density in the population
  Dens <- n.groups/((Xu - Xl) * (Yu - Yl))
  
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
  p <- p0 * exp(-alpha1 * D * D)
  
  ##-- Set-up observation matrix y:
  ##-- the dimensions are n.traps*K and it contains the number of wolves 
  ##-- photographed together for each trap and occasion
  V <- array(NA, c(n.groups+1,n.traps,n.occasions))
  dimnames(V) <- list("groups" = c(1:n.groups,"undetected"),
                      "traps" = paste("trap", 1:n.traps),
                      "occasions" = 1:n.occasions)
  Y <- matrix(0, nrow = n.traps, ncol = n.occasions)
  dimnames(Y) <- list("traps" = paste("trap", 1:n.traps),
                      "occasions" = 1:n.occasions)
  ##-- Loop over occasions
  for(k in 1:n.occasions){
    
    ##-- Sample group visit (or absence of visit)
    V[ , ,k] <- apply(p, 2, function(thisCol){
      sum.p <- sum(thisCol)
      risk <- 1 - exp(-sum.p)
      thisCol <- risk * thisCol/sum.p
      newP <- c(thisCol, 1 - risk) 
      out <- rmulti(n = 1,
                    size = 1,
                    prob = newP)
      return(out)
    })
    
    ##-- Sample number of wolves detected together for each visit
    for(j in 1:n.traps){
      ##-- Identify which group visited detector j (if any)
      whichGroup <- which(V[1:n.groups,j,k] > 0)
      ##-- Sample number of wolves detected together
      ##-- (at least one wolf detected)
      if(length(whichGroup)>0){
        Y[j,k] <-  extraDistr::rtbinom(n = 1,
                                       size = GS[whichGroup],
                                       prob = alpha,
                                       a = 0)
      }#if
    }#j
  }#k
  
  ##-- Output
  return(list( "Y" = Y,
               "n.groups" = n.groups,
               "N" = sum(GS),
               "n.traps" = n.traps,
               "n.occasions" = n.occasions,
               "groupSize" = GS, 
               "S" = S,
               "trapCoords" = trapCoords,
               "habitatBoundaries" = list("x" = c(Xl,Xu),
                                          "y" = c(Yl,Yu)),
               "visits" = V[1:n.groups, , ]))
}


testSim <- sim.SGM( n.groups = 15,         ## Number of groups
                    n.occasions = 30,      ## Number of detection occasions
                    p0 = 0.1,              ## Baseline probability of visit 
                    sigma = 0.8,           ## Scale parameter
                    lambda = 7,            ## Mean group size in the population
                    alpha = 0.85,          ## Cohesion of the group (probability to detect X ids together)
                    traps_dim = c(8,10))


testSim$N
##-- Distribution of number of pictures per trap
hist(apply(testSim$visits,2,sum),breaks = seq(-0.5,10.5,1))
##-- Distribution of numbers of pictures per group
hist(apply(testSim$visits,1,sum),breaks = seq(-0.5,20.5,1))



## ----------------------------------------------------------------------------
## ----- 2. NIMBLE MODEL ------
## -----   2.1. CUSTOM NIMBLE FUNCTIONS -----
## This function calculates the probability of detection of a group using a half-normal detection function
## It uses an indicator to make calculation faster (skips calculations for augmented groups)
calculateP <- nimbleFunction(
  run = function( p0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  indicator = double(0)){
    ## Specify return type
    returnType(double(1))
    
    ## Identify dimensions
    numTraps <- nimDim(trapCoords)[1]
    
    ## Calculate trap-specific detection rates
    p <- nimNumeric(length = numTraps)
    
    ## Shortcut if individual is not available for detection
    if (indicator < 1){return(p)}
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Loop over traps
    for(j in 1:numTraps){
      d2 <- pow(trapCoords[j,1] - s[1], 2) + pow(trapCoords[j,2] - s[2], 2)
      p[j] <- p0 * exp(alpha * d2)
    }#j
    return(p)
  })


## This custom distribution is an attempt to represent the detection process of
## wolf packs at a given camera trap.
## 1 - if x = 0; it returns the probability to detect 0 groups given the group-specific probabilities
## 2-  if x > 0; it returns the probability to detect x individuals out of a group of size[g] for the different groups
## using a truncated binomial x ~ T(dbinom( ))
dgroupDet <- nimbleFunction(
  run = function( x = double(0),           ## data: Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  p = double(1),           ## Group-specific detection probabilities
                  alpha = double(0),       ## Cohesion (i.e. prob of detection given visit)
                  indicator = double(1),   ## Group-specific augmentation indicator
                  log = integer(0, default = 0)
  ){
    ##-- Return type
    returnType(double(0))
    
    ##-- Derive probability of 0 visit:
    numGroups <- length(p)
    sumP <- sum(p[1:numGroups])
    pNULL<- exp(-sumP)
    risk <- 1 - pNULL
    
    ##-- if n = 0, return probability of 0 visit :
    if(x == 0){
      if(log == 0){return(pNULL)}else{return(log(pNULL))}
    }
    
    ##-- If x > 0; loop over groups and calculate the probability
    ##-- to detect X individuals simultaneously given the group size 
    ##-- and group-specific probability of visit :
    totalProb <- 0
    for(g in 1:numGroups){
      if(indicator[g] > 0){
        if(x <= size[g]){
          ##-- Probability of visit by group[g]
          probThisGroup <- risk * p[g]/sumP
          ##-- Probability of x wolves from group[g] detected together 
          logNormCst <- pbinom(q = 0,
                               size = size[g],
                               prob = alpha,
                               lower.tail = 0,
                               log.p = 1)
          
          logProbThisNum <- dbinom(x = x,
                                   size = size[g],
                                   prob = alpha,
                                   log = 1) -  logNormCst
          
          ##-- Probability of the data
          totalProb <- totalProb + probThisGroup * exp(logProbThisNum)
        }#if
      }
    }#g
    if(log) return(log(totalProb)) else return(totalProb)
  })




## -----   2.2. SPATIAL GROUP MODEL DEFINITION  -----
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  for(g in 1:G){
    s[g,1] ~ dunif(xlim[1], xlim[2])
    s[g,2] ~ dunif(ylim[1], ylim[2])
  }#g
  
  
  ##---- DEMOGRAPHIC PROCESS 
  psi ~ dunif(0,1)
  lambda ~ dunif(0,10)
  
  for(g in 1:G){
    z[g] ~ dbern(psi)
    groupSize[g] ~ T(dpois(lambda),1,20 )
  }#g
  
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  alpha ~ dunif(0,1) ## cohesion parameter (proportion of the group detected together)
  
  ##-- Calculate individual detection prob. at all traps
  for(i in 1:G){ 
    p[i,1:J] <- calculateP(
      p0 = p0,
      sigma = sigma,
      s = s[i,1:2],
      trapCoords = trapCoords[1:J,1:2],
      indicator = z[i])
  }#i
  
  ##-- Detection probability of one group (= probability of visit for now)
  for(j in 1:J){
    for(k in 1:K){
      y[j,k] ~ dgroupDet( 
        size = groupSize[1:G],        
        p = p[1:G,j],
        alpha = alpha,
        indicator = z[1:G])
    }#k
  }#j
  
  
  ##-- DERIVED PARAMETERS 
  n.groups <- sum(z[1:G])
  N <- sum(z[1:G]*groupSize[1:G])
})



## -----   2.3. BUNDLE DATA ------
##-- Nimble data
nimData <- list( y = testSim$Y)

##-- Nimble constants
nimConstants <- list( G = testSim$n.groups*2,
                      J = testSim$n.traps,
                      K = testSim$n.occasions,
                      trapCoords = as.matrix(testSim$trapCoords),
                      xlim = testSim$habitatBoundaries$x,
                      ylim = testSim$habitatBoundaries$y)


##-- Nimble inits
s.init <- rbind(testSim$S,
                testSim$S)
z.init <- c(rep(1,testSim$n.groups),
            rep(0,testSim$n.groups))
groupSize.init <- c(testSim$groupSize,
                    testSim$groupSize)
nimInits <- list( p0 = runif(1,0,1),
                  p = matrix(0.05,
                             nrow = nimConstants$G,
                             ncol = nimConstants$J),
                  psi = 0.5,
                  sigma = 0.5,           
                  lambda = 3,            
                  alpha = 0.75,
                  groupSize = groupSize.init,
                  z = z.init,
                  s = s.init)



## -----   2.4. FIT MODEL ------
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
Cmodel$calculate()
modelConf <- configureMCMC( model,
                            monitors = c( "N", "n.groups", "lambda",
                                               "psi", "p0", "sigma", "alpha"))
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
system.time(nimOutput <- runMCMC(CmodelMCMC,
                            niter = 1000,
                            nchains = 1,
                            nburnin = 0,
                            samplesAsCodaMCMC = T))
plot(nimOutput)



## ----------------------------------------------------------------------------