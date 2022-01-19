################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[.]_z[psi,rho]_y[p0,sigma] ---------------------- #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(rgdal)
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


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(file.path(getwd(),"Source"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "sim.0.0"

## HABITAT SPECIFICATIONS
habitat = list( country =  c("SWE","NOR"),
                resolution = 10000, 
                buffer = 60000)

## NGS DATA SPECIFICATIONS
dna = NA #list( sex = c("female","male")) 

## DETECTORS SPECIFICATIONS
detectors = list( detSubResolution = 1000,
                  detResolution = 10000)

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(simDir))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(simDir, modelName))){dir.create(file.path(simDir, modelName))}



## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
## Study area grid
studyAreaGrid <- readOGR(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- spTransform(x = studyAreaGrid, CRSobj = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
studyArea <- aggregate(studyAreaGrid)

## NEED COUNTRY LEVEL POLYGON TO CUT OUT SEA AND OTHER NON-HABITAT
countries <- readOGR(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## Plot check
plot(studyAreaGrid, col = "gray80")
plot(countries, col = "gray80", add = T)
plot(studyAreaGrid, col = adjustcolor("forestgreen",0.2), add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 1, ], col = "gray40", add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 2, ], col = "gray60", add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 3, ], col = "firebrick3", add = T)



## ------   2. SEARCH EFFORT DATA ------
## Load GPS search transects
transects <- readOGR(file.path(dataDir, "/GISData/Transects_Wolfalps20202021/wolfalps_transects_20202021.shp"))

## Plot check
plot(transects, col = "blue", add = T)



## ------   3. DNA DATA ------
# ngs <- read.csv(file.path(dataDir,
#                           "RovbaseData/ROVBASE DOWNLOAD 20191008/dna_wolf.csv"))
# ## Format dates
# ngs$Date <- as.POSIXct(strptime(ngs$Date, "%d.%m.%Y"))
# ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
# ngs$Month <- as.numeric(format(ngs$Date,"%m"))
# 
# ## Remove samples without coordinates
# ngs <- ngs[!is.na(ngs$North), ]
# 
# ## Remove samples without ID
# ngs <- ngs[!is.na(ngs$ID), ]
#
# ## Plot check
# plot(ngs, add = T)



## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. HABITAT ------
## ------     1.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (study area + buffer)  
habitat$polygon <- rgeos::gBuffer( studyArea,
                                   width = habitat$buffer)
## Remove un-suitable habitat (e.g. seas)
habitat$polygon <- rgeos::gIntersection( aggregate(countries),
                                         habitat$polygon,
                                         drop_lower_td = TRUE) 

##---- Create Habitat raster (keep raster cells with > 50% habitat) 
habitat$raster <- raster(extent(habitat$polygon))

##---- Trick to get % of habitat
res(habitat$raster) <- habitat$resolution/10
habitat$raster <- fasterize( st_as_sf(habitat$polygon),
                             habitat$raster)
habitat$raster[is.na(habitat$raster)]<- 0
habitat$raster <- aggregate(habitat$raster, fact = 10, sum, na.rm = TRUE)
habitat$raster[habitat$raster < 50] <- 0
habitat$raster[habitat$raster >= 50] <- 1

##---- Create Habitat matrix of cell ID 
habitat$grid <- habitat$raster
## Identify suitable habitat cells
isHab <- which(habitat$raster[]==1)
## Cell ID starts from the top left corner, increments left to right and up to down
habitat$grid[isHab] <- 1:length(isHab)
## Convert to matrix
habitat$grid <- as.matrix(habitat$grid)
habitat$matrix <- habitat$grid
habitat$matrix[habitat$matrix > 0] <- 1 

##---- Obtain xy coordinates of habitat cells   
habitat$coords <- coordinates(habitat$raster)[isHab, ]
dimnames(habitat$coords) <- list(1:length(isHab), c("x","y"))
# habitat.sp <- SpatialPointsDataFrame(data.frame(habCoords[ ,c("x","y")]),
#                                      data = data.frame(habCoords),
#                                      proj4string = CRS(projection(habitat)))
## Retrieve habitat windows boundaries
habitat$lowerCoords  <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)

##---- Visual plotting to check if everything is right  
plot(habitat$raster)					
plot(countries, add = T)
plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
plot(studyArea, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 1))



## ------     1.2. HABITAT COVARIATES ------
##---- Simulate 2 random habitat covariates
habitat$covariates <- cbind.data.frame(
  elev = rnorm(n = dim(habitat$coords)[1]),
  cover = rnorm(n = dim(habitat$coords)[1]))



## ------   2. DETECTORS ------
## ------     2.1. DETECTORS CHARACTERISTICS ------
##---- Create a raster of searched grid cells
detectors$raster <- raster(extent(transects),
                           resolution = detectors$detResolution,
                           crs = projection(transects))

##---- Extract transect length per detector
transects.p <- SampleTrack( x = transects,
                            d = 50)
detectors$raster <- rasterize( transects.p,
                               detectors$raster,
                               field = "LID",
                               fun = "count")
detectors$raster[] <- detectors$raster[]*20
isDet <- which(!is.na(detectors$raster[])) 
detectors$coords <- cbind.data.frame( "x" = coordinates(detectors$raster)[isDet,1],
                                      "y" = coordinates(detectors$raster)[isDet,2],
                                      "length" = detectors$raster[isDet])
detectors$n.detectors <- dim(detectors$coords)[1] 

##---- Visual plotting to check if everything is right  
plot(detectors$raster, add = T)
plot(transects, col = "blue", add = T)



## ------     2.2. DETECTOR COVARIATES ------





## ------   3. RESCALE HABITAT & DETECTORS ------
scaledCoords <- scaleCoordsToHabitatGrid(
  coordsData = detectors$coords,
  coordsHabitatGridCenter = habitat$coords)

detectors$scaledCoords <- scaledCoords$coordsDataScaled
habitat$scaledCoords <- scaledCoords$coordsHabitatGridCenterScaled
habitat$loScaledCoords <- habitat$scaledCoords - 0.5
habitat$upScaledCoords <- habitat$scaledCoords + 0.5

localObjects <- getLocalObjects(
  habitatMask = habitat$matrix,
  coords = detectors$scaledCoords[ ,1:2],
  dmax = 15)



## ------   4. DETECTION DATA ------
## ------     4.1. DETECTION MATRIX : y ------
# y.ar <- MakeY( myData = myData.alive$myData.sp,
#                myDetectors = myDetectors$main.detector.sp,
#                method = "Bernoulli",
#                returnIdvector = TRUE)

## SIMULATE DATA FOR NOW?

## ------     4.2. INDIVIDUAL COVARIATES ------
## Sex
## Pack membership
## ...



## ------   5. DATA AUGMENTATION ------
# y <- MakeAugmentation( y = y,
#                        aug.factor = augFactor,
#                        replace.value = 0)
# 
# sex <- MakeAugmentation( y = sex,
#                          aug.factor = augFactor,
#                          replace.value = NA)


## -----------------------------------------------------------------------------
## ------ III. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##-- SPATIAL PROCESS  
  # for(c in 1:n.habCovs){betaHab[c] ~ dnorm(0.0,0.01)}
  # for(h in 1:n.habWindows){
  #   habIntensity[h] <- exp(inprod(betaHab[1:n.habCovs],hab.covs[h,1:n.habCovs]))
  #   }
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##-- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)   
  rho ~ dunif(0,1)
  
  for(i in 1:n.individuals){ 
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(rho)
  }#i 								
  
  
  ##-- DETECTION PROCESS 
  p0[1] ~ dunif(0,1)
  p0[2] ~ dunif(0,1)
  sigma[1] ~ dunif(0,6)
  sigma[2] ~ dunif(0,6)
  
  # betaDet ~ dnorm(0.0,0.01)
  # logit(p0Traps[1,1:n.detectors]) <- logit(p0[1]) + betaDet * searchEffort[1:n.detectors]
  # logit(p0Traps[2,1:n.detectors]) <- logit(p0[2]) + betaDet * det.covs[1:n.detectors]
  
  for(i in 1:n.individuals){
    y[i,1:n.maxDets] ~ dbinomLocal_normal(
      size = size[1:n.detectors],
      p0 = p0[sex[i]+1],
      sigma = sigma[sex[i]+1],
      s = sxy[i,1:2],
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


## ------   2. DATA SIMULATION ------
simConstants <- list( n.individuals = 500,
                      n.habWindows = habitat$n.HabWindows,
                      # n.habCovs = 2,
                      n.detectors = detectors$n.detectors, 
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      n.maxDets = 40,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

simData <- list( lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords, 
                 habIntensity = rep(1,habitat$n.HabWindows),
                 size = rep(1,detectors$n.detectors),
                 detCoords = detectors$scaledCoords[,1:2],
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid = localObjects$habitatGrid)

simInits <- list( "psi" = 0.8,
                  "rho" = 0.5,
                  "p0" = c(0.06,0.04),
                  "sigma" = c(0.8,1)) 

##---- Create the nimble model object
simModel <- nimbleModel( code = modelCode,
                         constants = simConstants,
                         inits = simInits,
                         data = simData,
                         check = FALSE,
                         calculate = FALSE) 

##---- Get nodes to simulate
nodesToSim <- simModel$getDependencies(c(names(simInits),"habIntensity"),
                                       self = F,
                                       includeData = F,
                                       downstream = T,
                                       returnScalarComponents = T)

nodesToSim.simple <- unique(sapply(strsplit(nodesToSim, "\\["), "[", 1))
nodesToSim.simple

##---- Simulate from NIMBLE 
simModel$simulate(nodesToSim)
simModel$calculate()



## ------   3. BUNDLE NIMBLE DATA/INITIAL VALUES/... ------
nimConstants <- simConstants 
nimData <- simData
nimInits <- simInits

##---- Identify which individuals were detected
whichDet <- simModel$y[ ,1] > 0
whichNotDet <- simModel$y[ ,1] == 0

##---- Add simulated data 
nimData$y <- simModel$y

z <- z.init <- simModel$z
z[whichNotDet] <- NA
z.init[whichDet] <- NA
nimData$z <- z
nimInits$z <- z.init

sex <- sex.init <- simModel$sex
sex[whichNotDet] <- NA
sex.init[whichDet] <- NA
nimData$sex <- sex
nimInits$sex <- sex.init

nimInits$sxy <- simModel$sxy

nimParams <- c("N", 
               "p0", "sigma", "psi", "rho")
               #"z", "sxy")



## ------   4. FIT MODEL -----
##---- Set-up NIMBLE model
nimModel <- nimbleModel( code = modelCode,
                         constants = nimConstants,
                         data = nimData,
                         inits = nimInits,
                         check = F,       
                         calculate = F)  
conf <- configureMCMC( model = nimModel,
                       monitors = nimParams,
                       thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
for(c in 1:2){
  runMCMCbites( mcmc = Cmcmc,
                bite.size = 500,
                bite.number = 5,
                path = file.path(simDir, modelName, paste0("output/chain",c)))  
}#c

##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(simDir, modelName, "output"),
                               burnin = 2)

##---- Save MCMC output & plot
save( nimOutput,
      file = file.path(simDir, modelName, "mcmc_output.RData"))

pdf(file = file.path(simDir, modelName, "traceplots.pdf"))
plot(nimOutput)
graphics.off()


## -----------------------------------------------------------------------------
