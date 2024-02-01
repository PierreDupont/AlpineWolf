# I. DATA SET UP 
## Create habitat grid cell coordinates
xBounds <- c(0.5,25.5)
yBounds <- c(0.5,50.5)
habResolution <- 1
xCoords <- seq(xBounds[1],xBounds[2],habResolution)
yCoords <- seq(yBounds[1],yBounds[2],habResolution)
coordsHabitatGridCenter <- expand.grid(
  list(x = xCoords,
       y = yCoords))
colnames(coordsHabitatGridCenter) <- c("x","y")
numHabWindows <- nrow(coordsHabitatGridCenter)

## Create trap coordinates
bufferWidth <- 2
trapResolution <- 0.7
trapCoords <- expand.grid(
  list(x = seq(xBounds[1]+bufferWidth,xBounds[2]-bufferWidth,trapResolution),
       y = seq(yBounds[1]+bufferWidth,yBounds[2]-bufferWidth,trapResolution)))
colnames(trapCoords) <- c("x","y")
numTraps <- nrow(trapCoords)
oper <- rep(1,numTraps)

## RESCALE COORDINATES 
scaledCoords <- scaleCoordsToHabitatGrid( coordsData = trapCoords,
                                          coordsHabitatGridCenter = coordsHabitatGridCenter)
scaledtrapCoords <- scaledCoords$coordsDataScaled
habitatMask <- matrix(1, nrow = length(yCoords), ncol = length(xCoords), byrow = TRUE)


## CREATE LOCAL OBJECTS 
localTraps <- getLocalObjects( habitatMask = habitatMask,
                               coords = scaledtrapCoords,
                               dmax = 10,
                               resizeFactor = 1,
                               plot.check = TRUE)

## simulate population
N <- 50
M <- 150
indicator <- c(rep(1,N),rep(0,M-N))
s <- cbind(runif(n = M, min = xBounds[1], max = xBounds[2]),
           runif(n = M, min = yBounds[1], max = yBounds[2]))
colnames(s) <- c("x","y")
scaledS <- scaleCoordsToHabitatGrid( coordsData = s,
                                    coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled

## PLOT CHECK (RESCALE SCALE)
plot(scaledCoords$coordsHabitatGridCenterScaled[,"y"]~scaledCoords$coordsHabitatGridCenterScaled[,"x"]) 
points(scaledtrapCoords[,"y"]~scaledtrapCoords[,"x"], col = "blue",pch=3) 
points(scaledS[,"y"]~scaledS[,"x"],col=adjustcolor("red",0.3),pch=16)
points(scaledS[indicator>0,"y"]~scaledS[indicator>0,"x"],col="red",pch=16) 


## PARAMETERS
lambda0 <- 1
sigma <- 2

# III. USING THE RANDOM GENERATION FUNCTION 
y <- rpoisLocalSC_normal( n = 1,
                          lambda0 = lambda0,
                          sigma = sigma,
                          s = scaledS,
                          trapCoords = scaledtrapCoords,
                          oper = oper,
                          localTrapsIndices = localTraps$localIndices,
                          localTrapsNum = localTraps$numLocalIndices,
                          resizeFactor = 1,
                          habitatGrid = localTraps$habitatGrid,
                          indicator = indicator)

## Plot number of detections per trap
points( scaledtrapCoords[,"y"]~scaledtrapCoords[,"x"],
        pch = 21,
        bg = adjustcolor("black", 0.5),
        col = "black",
        cex = y/2)


# II. USING THE DENSITY FUNCTION 
# OPTION 1: USING THE RANDOM GENERATION FUNCTIONNALITY 
# system.time(
#   dpoisSC_normal2( x = y,
#                    lambda0 = lambda0,
#                    sigma = sigma,
#                    s = s,
#                    trapCoords = trapCoords,
#                    oper = oper,
#                    indicator = indicator,
#                    log = T)
# )


system.time(
  dpoisSC_normal( x = y,
                  lambda0 = lambda0,
                  sigma = sigma,
                  s = s,
                  trapCoords = trapCoords,
                  oper = oper,
                  indicator = indicator,
                  log = T)
)

system.time(
  dpoisLocalSC_normal( x = y,
                       lambda0 = lambda0,
                       sigma = sigma,
                       s = scaledS,
                       trapCoords = scaledtrapCoords,
                       oper = oper,
                       localTrapsIndices = localTraps$localIndices,
                       localTrapsNum = localTraps$numLocalIndices,
                       resizeFactor = 1,
                       habitatGrid = localTraps$habitatGrid,
                       indicator = indicator,
                       log = T)
)



###-----------------------------------------------------------------------------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  # psiRJ ~ dunif(0, 1) # inclusion prob
  # for(c in 1:n.habCovs){
  #   betaHab.raw[c] ~ dnorm(0.0,0.01)
  #   zRJ[c] ~ dbern(psiRJ)
  #   betaHab[c] <- betaHab.raw[c] * zRJ[c]
  # }#c
  
  ##-- Intensity of the AC distribution point process
  # habIntensity[1:n.habWindows] <- exp(
  #   hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ##-- Loop over individuals
  for(i in 1:M) {
    s[i,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##---- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)

  for(i in 1:M){ 
    z[i] ~ dbern(psi)
  }#i 								
  
  
  ##---- DETECTION PROCESS 
  ## Sample number of days active for cameras with unknown values
  lambda.oper ~ dgamma(1,1)
  for(j in 1:J){
    oper[j] ~ T(dpois(lambda.oper),1, )
  }#j
  
  ## Detection parameters priors
  ## Use mixture of normal distributions based on SCR results instead?
  ## Use proportion of each sex/status as weights for the mixture
  ## Use mean and sd from posteriors as mu and sigma for Normal distributions
  sigma ~ dunif(0,10)
  lambda0 ~ dunif(0,10)
  
  n[1:J] ~ dpoisLocalSC_normal(
      lambda0 = lambda0,
      sigma = sigma,
      s = s[1:M,1:2],
      trapCoords = trapCoords[1:J,1:2],
      oper = oper[1:J],
      localTrapsIndices = localTrapsIndices[1:n.habWindows,1:n.localIndicesMax],
      localTrapsNum = localTrapsNum[1:n.habWindows],
      resizeFactor = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      indicator = z[1:M])

  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
  D <- N/area
})


  
  

