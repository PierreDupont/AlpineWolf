
## The goal is to further develop the spatial count model (AKA the unmarked SCR) 
## based on camera trap data for group-living species. The general idea is to 
## leverage the information contained in camera-trap pictures when multiple 
## individuals are detected together.

## The assumptions are:
## - all individuals detected together belong to the same group.
## - all individuals in a group share the same AC and home-range size.


## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list = ls())

## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")

modelName = "simulation_SG.0.1"
thisDir <- file.path(simDir, "SpatialGroup", modelName)

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir, recursive = T)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}



## ------ IMPORT REQUIRED LIBRARIES & FUNCTIONS ------
library(nimble)
library(R.utils)
library(coda)
library(ggplot2)
library(hrbrthemes)

# library(nimbleSCR)
sourceDirectory(file.path(gitDir,"Source"), modifiedOnly=FALSE)


## -----  CUSTOM DATA SIMULATION FUNCTION -----

## Adapted sim function from Catherine Sun
## changed to simulate group activity centers and group detections instead of individuals
sim.SGM <- function( n.groups = 10,       ## Number of groups in the population
                     n.occasions = 20,    ## Number of occasions
                     p0 = 0.75,           ## Baseline probability of visit 
                     sigma = 0.5,         ## Scale parameter
                     lambda = 3,          ## Mean group size in the population
                     alpha = 0.75,        ## Cohesion of the group (probability to detect X ids together)
                     traps_dim = c(5,5),  ## Detector grid size
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




## -----------------------------------------------------------------------------
## ------ I. GENERATE THE PARAMETER SPACE TO BE TESTED -----
param.space <- MakeParameterSpace( 
  list.param = list( "n.groups" = c(3,10,20),
                     "lambda" = c(1,5,15),
                     "n.occasions" = c(1,10,30),
                     "p0" = 0.2,
                     "sigma" = 1,
                     "alpha" = 0.85),
  n.rep = 20)



## -----------------------------------------------------------------------------
## ------ II. DATA SIMULATION ------
for(r in 1:dim(param.space)[1]){
  
  ### ------ 1. DATA SIMULATION -----
  
  print(r)
  
  paramspace <- param.space[r, ]
  
  testSim <- sim.SGM( n.groups = param.space$n.groups[r],       ## Number of groups
                      n.occasions = param.space$n.occasions[r], ## Number of detection occasions
                      p0 = param.space$p0[r],                   ## Baseline probability of visit 
                      sigma = param.space$sigma[r],             ## Scale parameter
                      lambda = param.space$lambda[r],           ## Mean group size in the population
                      alpha = param.space$alpha[r],             ## Cohesion of the group (probability to detect X ids together)
                      traps_dim = c(8,10))
  
  
  ## ---------------------------------------------------------------------------
  
  ## ----- 2. SPATIAL GROUP MODEL ------
  ##-- Nimble model
  modelCode <- nimbleCode({
    ##---- SPATIAL PROCESS  
    for(g in 1:G){
      s[g,1] ~ dunif(xlim[1], xlim[2])
      s[g,2] ~ dunif(ylim[1], ylim[2])
    }#g
    
    
    ##---- DEMOGRAPHIC PROCESS 
    psi ~ dunif(0,1)
    lambda ~ dgamma(1,1)
    
    for(g in 1:G){
      z[g] ~ dbern(psi)
      groupSize[g] ~ T(dpois(lambda),1,)
    }#g
    
    
    ##---- DETECTION PROCESS
    sigma ~ dunif(0,10)
    p0 ~ dunif(0,1)
    alpha ~ dunif(0,1) ## cohesion parameter (proportion of the group detected together)
    
    ##-- Calculate individual detection prob. at all traps
    for(g in 1:G){ 
      p[g,1:J] <- calculateP(
        p0 = p0,
        sigma = sigma,
        s = s[g,1:2],
        trapCoords = trapCoords[1:J,1:2],
        indicator = z[g])
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
  
  ##-- Nimble data
  nimData <- list( y = testSim$Y)
  
  ##-- Nimble constants
  nimConstants <- list( G = testSim$n.groups * 2,
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
  
  ##-- Nimble parameters
  nimParams <- c( "N", "n.groups", "lambda",
                  "psi", "p0", "sigma", "alpha")
  nimParams2 <- c("s","z")
  
  ##-- Save input
  thiFile <- file.path( thisDir, "input",
                        paste0( "SpatialGroup_", param.space$set_ID[r],
                                "Rep", param.space$rep_ID[r],
                                ".RData"))
  save( modelCode,
        nimData,
        nimConstants,
        nimInits,
        nimParams,
        nimParams2,
        paramspace,
        testSim,
        file = thiFile)
  
  
  
  
  ## ---------------------------------------------------------------------------
  
  ## ----- 3. SPATIAL COUNT MODEL ------
  ##-- Nimble model
  modelCode <- nimbleCode({
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
    lambda0 ~ dgamma(1,1)
    
    y[1:J] ~ dpoisSC_normal(
      lambda0 = lambda0,
      sigma = sigma,
      s = s[1:M,1:2],
      oper = oper[1:J],
      trapCoords = trapCoords[1:J,1:2],
      indicator = z[1:M])
    
    
    ##-- DERIVED PARAMETERS 
    N <- sum(z[1:M])
  })
  
  ##-- Nimble data
  nimData <- list( y = rowSums(testSim$Y))
  
  ##-- Nimble constants
  nimConstants <- list( M = testSim$N * 2,
                        J = testSim$n.traps,
                        oper = rep(testSim$n.occasions, testSim$n.traps),
                        trapCoords = as.matrix(testSim$trapCoords),
                        xlim = testSim$habitatBoundaries$x,
                        ylim = testSim$habitatBoundaries$y)
  
  ##-- Nimble inits
  s.init <- cbind( runif( n = nimConstants$M,
                          min = nimConstants$xlim[1], 
                          max = nimConstants$xlim[2]),
                   runif( n = nimConstants$M,
                          min = nimConstants$ylim[1], 
                          max = nimConstants$ylim[2]))
  
  z.init <- c(rep(1,testSim$N),
              rep(0,testSim$N))
  
  nimInits <- list( lambda0 = runif(1,0,10),
                    psi = 0.5,
                    sigma = 0.5,           
                    z = z.init,
                    s = s.init)
  
  ##-- Nimble parameters
  nimParams <- c("N","lambda0","psi","sigma")
  nimParams2 <- c("s","z")
  
  ##-- Save input
  thiFile <- file.path( thisDir, "input",
                        paste0( "SpatialCount_", param.space$set_ID[r],
                                "Rep", param.space$rep_ID[r],
                                ".RData"))
  save( modelCode,
        nimData,
        nimConstants,
        nimInits,
        nimParams,
        nimParams2,
        paramspace,
        testSim,
        file = thiFile)
}#r



## -----------------------------------------------------------------------------
## ------ III. MODEL FITTING ------

r <- sample(x = nrow(param.space), size = 1)

load( file.path( thisDir, "input",
                             paste0( "SpatialGroup_", param.space$set_ID[r],
                                     "Rep", param.space$rep_ID[r],
                                     ".RData")))
model <- nimbleModel( code = modelCode,
                         constants = nimConstants,
                         data = nimData,
                         inits = nimInits)
Cmodel <- compileNimble(model)
Cmodel$calculate()
modelConf <- configureMCMC( model, monitors = nimParams)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble( modelMCMC, project = model)
system.time(nimOutput <- runMCMC( CmodelMCMC,
                               niter = 100,
                               nchains = 2,
                               samplesAsCodaMCMC = T))
plot(nimOutput)



## -----------------------------------------------------------------------------
## ------ IV. PROCESSING MODEL OUTPUTS ------

##-- List all input files
inFiles <- list.files( path = file.path(thisDir, "input"))

##-- Create a dataframe to store the results
output.df <- do.call( rbind,
                      lapply(inFiles,
                             function(x){
                               load(file.path(thisDir, "input", x))
                               cbind( "model" = ifelse(length(grep(pattern = "SpatialGroup",x = x)) > 0,
                                                       "SpatialGroup",
                                                       "SpatialCount"),
                                      paramspace,
                                      "N.sim" = testSim$N,
                                      "lambda.sim" = mean(testSim$groupSize),
                                      "mean.N" = NA,
                                      "sd.N" = NA,
                                      "mean.n.groups" = NA,
                                      "sd.n.groups" = NA,
                                      "mean.lambda" = NA,
                                      "sd.lambda" = NA,
                                      "mean.psi" = NA,
                                      "sd.psi" = NA,
                                      "mean.p0" = NA,
                                      "sd.p0" = NA,
                                      "mean.sigma" = NA,
                                      "sd.sigma" = NA,
                                      "mean.alpha" = NA,
                                      "sd.alpha" = NA,
                                      "mean.lambda0" = NA,
                                      "sd.lambda0" = NA)
                             }))



##-- List all output files
outFiles <- list.files(path = file.path(thisDir, "output"))

##-- Loop over inFiles
for(i in 1:length(inFiles)){
  ##-- check if the corresponding output exists
  thisFile <- paste0("NimOutput_", inFiles[i])
  
  ##-- if the output exists, process it and store results
  if(thisFile %in% outFiles){
    load(file.path(thisDir, "output", thisFile))
    out <- ProcessCodaOutput(nimOutput$samples)

    output.df$mean.N[i] <- out$mean$N
    output.df$sd.N[i] <- out$sd$N
    output.df$mean.sigma[i] <- out$mean$sigma
    output.df$sd.sigma[i] <- out$sd$sigma
    
    # if(output.df$model == "SpatialGroup"){
    #   output.df$mean.N <- out$mean$N
    #   out$mean$sigma
    # } else {
    #   output.df$mean.N <- out$mean$N
    #   
    # }
  }
}#i
output.df$RB.N <- (output.df$N.sim-output.df$mean.N)/output.df$N.sim

##-- VL violin plots

ggplot(data = output.df,
       aes(x = as.factor(lambda) ,
           y = RB.N,
           fill = as.character(model))) +
  geom_hline(aes(yintercept=0),  linetype ="solid", color = "red", linewidth=1) +
  geom_violin(trim = FALSE, linewidth =0) +
  stat_summary(fun = mean, geom = "point", color = "white", size = 0.5,
               position = position_dodge(0.9)) +
  facet_grid(output.df$n.groups ~ output.df$n.occasions) +
  labs(title=" ", fill ="Model Comparison", y = "Relative Bias", x = "Mean number of individuals per group") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=12),
        axis.text.x = element_text(size=12,angle=0,hjust=1),
        # axis.title.x = element_text(size=12,angle=0,hjust=1),
        axis.text.y = element_text(size=12),
        # axis.title.y = element_text(size=12,angle=0,hjust=1),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position="bottom", legend.box = "horizontal") 

# ggsave("modelSCSG.png", width = 25, height = 15, units = "cm")









