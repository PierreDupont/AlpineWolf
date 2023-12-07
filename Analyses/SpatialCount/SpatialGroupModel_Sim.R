
library(extraDistr)
library(dplyr)
library(ggplot2)
library(nimble)


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
  ##-- the dimensions are ntraps*K and it contains the number of wolves 
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
      whichGroup <- which(V[1:G,j,k] > 0)
      ##-- Sample number of wolves detected together
      if(length(whichGroup)>0){
        Y[j,k] <- rtbinom(n = 1,
                          size = GS[whichGroup],
                          prob = alpha,
                          a = 0)
      }#if
    }#j
  }#k
  
  ##-- Output
  return(list( "Y" = Y,
               "n.groups" = n.groups,
               "n.traps" = ntraps,
               "n.occasions" = n.occasions,
               "groupSize" = GS, 
               "S" = S,
               "trapCoords" = trapCoords,
               "habitatBoundaries" = list("x" = c(Xl,Xu),
                                          "y" = c(Yl,Yu)),
               "visits" = V))
}




## -----   1.2. SIMULATE DATA ------
testSim <- sim.SGM()




## ----------------------------------------------------------------------------
## ----- 2. MODEL DEFINITION ------
## -----  2.1. CUSTOM FUNCTION ------
dbinom_vector_truncated <- nimbleFunction(
  run = function( x = double(0),           ## Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  prob = double(1),        ## Group-specific detection probs (includes visits)
                  log = integer(0, default = 0)
  ){
    returnType(double(0))
    
    ##-- Shortcut if the number of pictures taken is 0
    if(x == 0){
      if(sum(prob) == 0){
        if (log == 0) 
          return(1)
        else return(0)
      }
      else {
        if (log == 0) 
          return(0)
        else return(-Inf)
      }
    }
    
    ##-- If detNums > 0; loop over groups and calculate the probability
    ##-- to get X individuals detected together given the group size and number of visits
    numGroups <- length(size)
    logProb <- 0
    
    for(g in 1:numGroups){
      thisProb <- dbinom(x,
                           prob = prob[g],
                           size = size[g],
                           log = TRUE)T(1, )
      logProb <- logProb + thisProb
    }#g
    if(log) return(logProb) else return(exp(logProb))
  })


## -----  2.2. NIMBLE MODEL ------

##-- SPATIAL GROUP MODEL (SGM)
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
    groupSize[g] ~ dpois(lambda)T(1, )
  }#g
  n.groups <- sum(z[1:G])
  N <- sum(z[1:G]*groupSize[1:G])
  D <- N/area
  
  
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  alpha ~ dunif(0,1) ## cohesion parameter (proportion of the group detected together)

  ##-- Detection probability of one group (= probability of visit for now)
  for(j in 1:n.traps){
    for(g in 1:G){
      d2[g,j] <- (s[g,1] - trapCoords[j,1])^2 + (s[g,2] - trapCoords[j,2])^2
      ##-- Group-specific probability of visit:
      p[g,j] <- p0 * exp(-d2[g,j] / (2*sigma^2)) * z[g] 
    }#g
    
    ##-- Derive probability of 0 visit:
    sum.p[j] <- sum(p[1:G,j])
    risk[j] <- 1 - exp(-sum.p[j])
    newP[1:G,j] <- risk[j] * p[1:G,j]/sum.p[j]
    newP[G1,j] <- 1 - risk[j]
    
    for(k in 1:n.occasions){
      ##-- Sample group visit :
      v[1:G1,j,k] ~ dmulti(
        size = 1,
        prob = newP[1:G1,j])
      
      ##-- Expected number of ids detected per picture per group:
      n[j,k] ~ dbinom_vector_truncated(
        size = groupSize[1:G],
        prob = alpha * v[1:G,j,k])
    }#k
  }#j
  
})




## ----------------------------------------------------------------------------
## ----- 3. MODEL FITTING ------
## -----  3.1. BUNDLE DATA ------

##-- Nimble data
nimData <- list( n = testSim$Y)


##-- Nimble constants
nimConstants <- list( G = testSim$n.groups*2,
                      G1 = testSim$n.groups*2+1,
                      n.traps = testSim$n.traps,
                      n.occasions = testSim$n.occasions,
                      trapCoords = testSim$trapCoords,
                      xlim = testSim$habitatBoundaries$x,
                      ylim =  testSim$habitatBoundaries$y)

##-- Nimble inits
z.init <- c(rep(1,testSim$n.groups),
            rep(0,testSim$n.groups))

v.init <- array(0,
                c(nimConstants$G1, nimConstants$n.traps, nimConstants$n.occasions))
v.init[1:testSim$n.groups, , ] <- testSim$visits[1:testSim$n.groups, , ]
v.init[nimConstants$G+1, , ] <- testSim$visits[testSim$n.groups+1, , ]

nimInits <- list( p0 = runif(1,0,1),
                  psi = 0.5,
                  sigma = 0.5,           ## Scale parameter
                  lambda = 3,            ## Mean group size in the population
                  alpha = 0.75,
                  z = z.init,
                  v = v.init)

##-- Nimble parameters
params <- c("N", "groupSize", "D", "n.groups", "lambda",
            "psi", "p0", "sigma", "alpha" )




## -----  3.2. FIT MODEL ------
model <- nimbleModel( code = SGM_model,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
out1 <- runMCMC(CmodelMCMC, niter = nadapt)
