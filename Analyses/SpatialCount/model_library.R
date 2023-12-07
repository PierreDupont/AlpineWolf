##------- UNMARKED MODEL -----
SC_model <- nimbleCode({
  
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] <- betaHab.raw[c] * zRJ[c]
  }#c
  
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    s[i,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid2[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  
  ##---- DEMOGRAPHIC PROCESS 
  psi ~ dunif(0,1)

  for(g in 1:n.individuals) {
    z[i] ~ dbern(psi)
  }#i

  N <- sum(z[1:n.individuals])
  D <- N/area
  
  
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  lambda0 ~ dgamma(1,1)
  
  ##-- Expected number of pictures per detector
  for(i in 1:n.individuals){
    for(j in 1:n.detectors){
      d2[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
      oneMinusP[i,j] <- p0 * exp(-d2[i,j] / (2*sigma^2)) * z[i]
      #lambda[i,j] <- lambda0 * exp(-d2[i,j] / (2*sigma^2)) * z[i]
    }#j
  }#g
  
  ##-- Binary/count detections per detector 
  for(j in 1:n.detectors){
    bigP[j] <- 1-prod(oneMinusP[1:n.individuals,j])
    y[j] ~ dbinom( size = n.occasions[j],
                   p = bigP[j])
    
    # bigLambda[j] <- sum(lambda[1:n.individuals,j])
    # n[j] ~ dpois(lambda = bigLambda[j])
  }#j

  ##-- This could me made faster by:
  ## 1- calculate 'bigP' in a custom function to avoid calculating distance when z = 0
  ## 2- use local evaluation; maybe not super efficient (?)
  })

  
  
  
##------- SPATIAL GROUP MODEL -----
## The idea of the SGM is to further develop the spatial count model 
## (AKA the unmarked SCR), to account for hte fact that multiple individuals are 
## detected simultaneously in one camera-trap picture.
SGM_model <- nimbleCode({
  
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] <- betaHab.raw[c]
  }#c
  
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(g in 1:n.groups){
    s[g,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid2[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#g
  
  
  
  ##---- DEMOGRAPHIC PROCESS 
  psi ~ dunif(0,1)
  lambdaGroup ~ dunif(0,10)
  
  for(g in 1:n.groups){
    z[g] ~ dbern(psi)
    groupSize[g] ~ T(dpois(lambdaGroup),1, )
  }#g
  G <- sum(z[1:n.groups])
  N <- sum(z[1:n.groups]*groupSize[1:n.groups])
  D <- N/area
    
  
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  p0 ~ dunif(0,1)
  alpha ~ dunif(0,1) ## cohesion parameter (proportion of the group detected together)
  #lambda0 ~ dgamma(1,1)
  
  ##-- Detection probability of one group (= probability of visit for now)
  for(j in 1:n.camtraps){
    for(g in 1:n.groups){
      d2[g,j] <- (s[g,1] - X[j,1])^2 + (s[g,2] - X[j,2])^2
      ##-- Group-specific probability of visit:
      p[g,j] <- p0 * exp(-d2[g,j] / (2*sigma^2)) * z[g] 
    }#g
    
    ##-- Derive probability of 0 visit:
    sum.p[j] <- sum(p[1:n.groups,j])
    risk[j] <- 1 - exp(-sum.p[j])
    newP[j,1:n.groups] <- risk[j] * p[1:n.groups,j]/sum.p[j]
    newP[j,n.groups+1] <- 1 - risk[j]
    
    ##-- Sample group visit :
    v[1:(n.groups+1),j] ~ dmulti( size = 1,
                                  prob = newP[j,1:(n.groups+1)])
    
    ##-- Expected number of ids detected per picture per group:
    n[j] ~ dbinom_vector_truncated( size = groupSize[1:n.groups],
                                    prob = alpha * v[1:n.groups,j])
  }#j
  
  
})
  
  


## ------ SPATIAL CAPTURE RECAPTURE MODEL ------
SCR_model <- nimbleCode({
  ##---- SPATIAL PROCESS  
  psiRJ ~ dunif(0, 1) # inclusion prob
  
  for(c in 1:n.habCovs){
    betaHab.raw[c] ~ dnorm(0.0,0.01)
    zRJ[c] ~ dbern(psiRJ)
    betaHab[c] <- betaHab.raw[c] * zRJ[c]
  }#c
  
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    s[i,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid2[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##---- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)
  rho ~ dunif(0,1)
  
  for(ss in 1:2){
    theta[1:n.states,ss] ~ ddirch(alpha[1:n.states,ss])
  }#ss
  
  for(i in 1:n.individuals){ 
    sex[i] ~ dbern(rho)
    z[i] ~ dbern(psi)
    status[i] ~ dcat(theta[1:n.states,sex[i]+1])
  }#i 								
  
  
  ##---- DETECTION PROCESS 
  for(c in 1:n.detCovs){
    betaDet[c] ~ dnorm(0.0,0.01)
  }
  
  for(s in 1:n.states){
    for(ss in 1:2){
      p0[s,ss] ~ dunif(0,0.5)
      sigma[s,ss] ~ dunif(0,12)
      logit(p0Traps[s,ss,1:n.detectors]) <- logit(p0[s,ss]) + 
        det.covs[1:n.detectors,1:n.detCovs] %*% betaDet[1:n.detCovs]
    }#ss
  }#s
  
  for(i in 1:n.individuals){
    y[i,1:n.maxDets] ~ dbinomLocal_normal(
      size = size[1:n.detectors],
      p0Traps = p0Traps[status[i],sex[i]+1,1:n.detectors],
      sigma = sigma[status[i],sex[i]+1],
      s = s[i,1:2],
      trapCoords = detCoords[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.habWindows,1:n.localIndicesMax],
      localTrapsNum = localTrapsNum[1:n.habWindows],
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      indicator = z[i],
      lengthYCombined = n.maxDets)
  }#i
  
  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:n.individuals])
})


## ------ SPATIAL INTEGRATED MODEL ------

SIM_model <- nimbleCode({
  
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab.raw[c] ~ dnorm(0.0,0.01)
    zRJ[c] ~ dbern(psiRJ)
    betaHab[c] <- betaHab.raw[c] * zRJ[c]
  }#c
  
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(g in 1:n.groups){
    s[g,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid2[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##---- DEMOGRAPHIC PROCESS 
  psi ~ dunif(0,1)
  lambda.group ~ dunif(0,10)
  
  for(g in 1:n.groups) {
    z[g] ~ dbern(psi)
    x[g] ~ dpois(lambda.group)
  }#g
  n.groups <- sum(z[1:n.groups])
  N <- sum(z[1:n.groups]*x[1:n.groups])
  D <- N/area
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0,10)
  lam0 ~ dunif(0,10)
  p0 ~ dunif(0,1)
  
  ##-- Expected number of pictures per detector
  for(g in 1:n.groups){
    for(j in 1:n.detectors){
      d2[g,j] <- (s[g,1] - X[j,1])^2 + (s[g,2] - X[j,2])^2
      lam[g,j] <- lam0 * exp(-d2[i,j] / (2*sigma^2)) * z[g] 
    }#j
  }#g
  
  ##-- Total number of pictures per detector (per occasion)
  for(j in 1:n.detectors){
    for(k in 1:n.occasions){
      bigLambda[j,k] <- sum(lam[1:n.groups,j]) * oper[j,k]
      n.pics[j,k] ~ dpois(bigLambda[j,k])
    }#k
  }#j
  
  ##-- Expected number of wolves per picture
  for(j in 1:n.detectors){
    for(k in 1:n.occasions){
      for(e in 1:n.events[j,k]){
        n.wolves[j,k,e] ~ dcustom(p0)
      }
    }
  }
  
  
})
