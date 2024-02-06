################################################################################
##### ----------------------- ALPINE WOLF SGM ---------------------------- #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(raster)
library(coda)
library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
library(lubridate)
library(stars)
library(extraDistr)
library(ggplot2)
#source("SpatialCount")

## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "sim_SGM_0.0"

## HABITAT SPECIFICATIONS
habitat = list( country =  c("SWE","NOR"),
                resolution = 10000, 
                buffer = 60000)

## NGS DATA SPECIFICATIONS
dna = NA #list( sex = c("female","male")) 

## DETECTORS SPECIFICATIONS
detectors = list( detSubResolution = 1000,
                  detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(simDir))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(simDir, modelName))){dir.create(file.path(simDir, modelName))}




## -----------------------------------------------------------------------------
## ----- 1. DATA SIMULATION ------
## -----   1.1. SIMULATION FUNCTION ------

# adapted sim function from oSCR, from Cat Sun
# changed to also output activity centers
sim.SGM <- function(n.groups = 10,                ## Number of groups
                    n.occasions = 20,                ## Number of detection occasions
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
               "n.traps" = n.traps,
               "n.occasions" = n.occasions,
               "groupSize" = GS, 
               "S" = S,
               "trapCoords" = trapCoords,
               "habitatBoundaries" = list("x" = c(Xl,Xu),
                                          "y" = c(Yl,Yu)),
               "visits" = V))
}




## -----   1.2. SIMULATE DATA ------
testSim <- sim.SGM( n.groups = 10,                ## Number of groups
                    n.occasions = 20,                ## Number of detection occasions
                    p0 = 0.75,             ## Baseline probability of visit 
                    sigma = 0.5,           ## Scale parameter
                    lambda = 3,            ## Mean group size in the population
                    alpha = 0.75,          ## Cohesion of the group (probability to detect X ids together)
                    traps_dim = c(5,5))




## ----------------------------------------------------------------------------
## ----- 2. MODEL DEFINITION ------
## -----  2.1. SPATIAL COUNT MODEL -----
SC_model <- nimbleCode({
  ##---- SPATIAL PROCESS 
  for(g in 1:G){
    s[g,1] ~ dunif(xlim[1], xlim[2])
    s[g,2] ~ dunif(ylim[1], ylim[2])
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
  D <- N/area
})



## -----  2.2. SPATIAL COUNT MODEL -----
##-- The idea of the SGM is to further develop the spatial count model 
##-- (AKA the unmarked SCR), to account for the fact that:
##--   1. multiple individuals from the same group are detected simultaneously 
##--   2. only one group is detected in a given picture
SGM_model <- nimbleCode({
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
    groupSize[g] ~ T(dpois(lambda),1, )
  }#g
  
  
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  alpha ~ dunif(0,1) ## cohesion parameter (proportion of the group detected together)
  
  ##-- Detection probability of one group (= probability of visit for now)
  for(j in 1:n.traps){
    for(k in 1:n.occasions){
      ##-- Sample group visit :
      n[j,k] ~ dcustom( 
        size = groupSize[1:G],        
        p0 = p0,
        sigma = sigma,
        alpha = alpha,
        s = s[1:G,1:2],
        trapCoords = trapCoords[j,1:2],
        indicator = z[1:G])
    }#k
  }#j
  
  ##-- DERIVED PARAMETERS 
  n.groups <- sum(z[1:G])
  N <- sum(z[1:G]*groupSize[1:G])
  D <- N/area
})




## ----------------------------------------------------------------------------
## ----- 3. MODEL FITTING ------
## -----  3.1. SPATIAL COUNT MODEL -----
## -----    3.1.1. BUNDLE DATA ------

##-- Nimble data
nimData <- list( n = testSim$Y,
                 numGroups = testSim$n.groups*2)

##-- Nimble constants
nimConstants <- list( G = testSim$n.groups*2,
                      n.traps = testSim$n.traps,
                      n.occasions = testSim$n.occasions,
                      trapCoords = testSim$trapCoords,
                      xlim = testSim$habitatBoundaries$x,
                      ylim =  testSim$habitatBoundaries$y)

##-- Nimble inits
z.init <- c(rep(1,testSim$n.groups),
            rep(0,testSim$n.groups))

groupSize.init <- c(testSim$groupSize,
                    testSim$groupSize)

nimInits <- list( p0 = runif(1,0,1),
                  psi = 0.5,
                  sigma = 0.5,           ## Scale parameter
                  lambda = 3,            ## Mean group size in the population
                  alpha = 0.75,
                  groupSize = groupSize.init,
                  z = z.init)#,
#v = v.init)

##-- Nimble parameters
params <- c("N", "groupSize", "D", "n.groups", "lambda",
            "psi", "p0", "sigma", "alpha" )



## -----    3.1.2. FIT MODEL ------
model <- nimbleModel( code = SGM_model,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
Cmodel$calculate()
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
out1 <- runMCMC(CmodelMCMC, niter = nadapt)

## -----  3.". SPATIAL GORUP MODEL -----
## -----    3.".1. BUNDLE DATA ------

##-- Nimble data
nimData <- list( n = testSim$Y,
                 numGroups = testSim$n.groups*2)

##-- Nimble constants
nimConstants <- list( G = testSim$n.groups*2,
                      n.traps = testSim$n.traps,
                      n.occasions = testSim$n.occasions,
                      trapCoords = testSim$trapCoords,
                      xlim = testSim$habitatBoundaries$x,
                      ylim =  testSim$habitatBoundaries$y)

##-- Nimble inits
z.init <- c(rep(1,testSim$n.groups),
            rep(0,testSim$n.groups))

groupSize.init <- c(testSim$groupSize,
                    testSim$groupSize)

nimInits <- list( p0 = runif(1,0,1),
                  psi = 0.5,
                  sigma = 0.5,           ## Scale parameter
                  lambda = 3,            ## Mean group size in the population
                  alpha = 0.75,
                  groupSize = groupSize.init,
                  z = z.init)#,
#v = v.init)

##-- Nimble parameters
params <- c("N", "groupSize", "D", "n.groups", "lambda",
            "psi", "p0", "sigma", "alpha" )



## -----    3..2. FIT MODEL ------
model <- nimbleModel( code = SGM_model,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
Cmodel$calculate()
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
out1 <- runMCMC(CmodelMCMC, niter = nadapt)
