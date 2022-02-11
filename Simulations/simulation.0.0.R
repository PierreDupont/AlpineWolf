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
library(lubridate)
library(stars)



## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)



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
                  detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(simDir))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(simDir, modelName))){dir.create(file.path(simDir, modelName))}




## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
## NEED COUNTRY LEVEL POLYGON TO CUT OUT SEA AND OTHER NON-HABITAT
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- st_union(studyAreaGrid)
# studyAreaGrid <- readOGR(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
# studyAreaGrid <- spTransform(x = studyAreaGrid, CRSobj = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
# studyArea <- aggregate(studyAreaGrid)

## Plot check
plot(studyArea, col = "gray80")
plot(countries["CNTR_CODE"], col = "gray80", add = T)
plot(studyAreaGrid["SCR"], add = T)




## ------   2. SEARCH EFFORT DATA ------
## Load GPS search transects
#transects <- read_sf(file.path(dataDir,"/GISData/Transects_Wolfalps20202021/wolfalps_transects_20202021.shp"))
#transects <- readOGR(file.path(dataDir, "/GISData/Transects_Wolfalps20202021/wolfalps_transects_20202021.shp"))
transects <- read_sf(file.path(dataDir,"/GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))

## Plot check
plot(transects, col = "red", add = T)



## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. HABITAT ------
##---- Create Habitat polygon (study area + buffer)  
# habitat$polygon <- rgeos::gBuffer( studyArea, width = habitat$buffer)
habitat$polygon <- st_buffer( studyArea, dist = habitat$buffer)

## Remove un-suitable habitat (e.g. seas)
habitat$polygon <- st_intersection( st_union(countries),
                                         habitat$polygon,
                                         drop_lower_td = TRUE) 

##---- Create Habitat raster (keep raster cells with > 50% habitat) 
# habitat$raster <- raster(extent(habitat$polygon))
habitat$raster <- raster(extent(st_bbox(habitat$polygon)))
## Trick to get % of habitat
res(habitat$raster) <- habitat$resolution/10
habitat$raster <- fasterize( st_as_sf(habitat$polygon),
                             habitat$raster)
habitat$raster[is.na(habitat$raster)]<- 0
habitat$raster <- aggregate(habitat$raster, fact = 10, sum, na.rm = TRUE)
habitat$raster[habitat$raster < 50] <- 0
habitat$raster[habitat$raster >= 50] <- 1

habitat$grid <- st_as_sf(rasterToPolygons(habitat$raster,
                                          fun = function(x)x == 1))
habitat$grid$id <- 1:nrow(habitat$grid)
habitat$grid <- habitat$grid[,c(2,3)]
st_crs(habitat$grid) <- st_crs(habitat$polygon)

##---- Create Habitat matrix of cell ID 
habitat$matrix <- habitat$raster
## Identify suitable habitat cells
isHab <- which(habitat$raster[]==1)
## Cell ID starts from the top left corner, increments left to right and up to down
habitat$matrix[isHab] <- 1:length(isHab)
## Convert to matrix
habitat$matrix <- as.matrix(habitat$matrix)
habitat$binary <- habitat$matrix
habitat$binary[habitat$binary > 0] <- 1 


##---- Obtain xy coordinates of habitat cells   
habitat$coords <- coordinates(habitat$raster)[isHab, ]
dimnames(habitat$coords) <- list(1:length(isHab), c("x","y"))
habitat$sp <- SpatialPoints(coords = habitat$coords,
                            proj4string = crs(habitat$polygon)) 


##---- Retrieve habitat windows corners
habitat$lowerCoords  <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)


##---- Visual plotting to check if everything is right  
plot(habitat$raster)					
plot(countries["CNTR_CODE"], col = NA, add = T)
plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
plot(studyArea, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 1))





## ------   2. DETECTORS ------
##---- Extract date and month for each transect
transects$Date <- as.numeric(str_sub(transects$Name, start= -6))
transects$Month <- as.numeric(str_sub(transects$Name, start= -4, end = -3))

##---- Create a detector grid from transects with the desired resolution
grid <- st_as_sf(st_make_grid(transects,detectors$detResolution))
grid$id <- 1:nrow(grid)
  
##---- Extract length and number of transects in each grid cell
intersection <- st_intersection(grid, transects) %>%
    mutate(LEN = st_length(.))  %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
              transect_N = length(unique(Date)))   ## Get total number of visits in each detector grid cell

##---- Store in detector grid
grid <- grid %>% 
  left_join(intersection, by = "id") %>% 
  filter(!is.na(transect_L))
plot(grid)
detectors$grid <- grid
  
##---- Extract detector coordinates
detectors$coords <- st_coordinates(st_centroid(detectors$grid))[ ,c("X","Y")]
dimnames(detectors$coords) <- list(1:nrow(detectors$grid),
                                   c("x","y"))

detectors$trials <- detectors$grid$transect_N

##---- Extract total number of detectors
detectors$n.detectors <- nrow(detectors$grid)




## ------   3. RESCALE HABITAT & DETECTORS ------
##---- Rescale coordinates
scaledCoords <- scaleCoordsToHabitatGrid(
  coordsData = detectors$coords,
  coordsHabitatGridCenter = habitat$coords)

detectors$scaledCoords <- scaledCoords$coordsDataScaled
habitat$scaledCoords <- scaledCoords$coordsHabitatGridCenterScaled
habitat$loScaledCoords <- habitat$scaledCoords - 0.5
habitat$upScaledCoords <- habitat$scaledCoords + 0.5


##---- Identify local detectors
localObjects <- getLocalObjects(
  habitatMask = habitat$binary,
  coords = scaledCoords$coordsDataScaled[ ,1:2],
  dmax = 21)



## -----------------------------------------------------------------------------
## ------ III. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i,1:2] ~ dbernppAC(
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
  }
  
  for(i in 1:n.individuals){ 
    sex[i] ~ dbern(rho)
    z[i] ~ dbern(psi)
    state[i] ~ dcat(theta[1:n.states,sex[i]+1])
  }#i 								
  
  
  ##---- DETECTION PROCESS 
  for(s in 1:n.states){
    for(ss in 1:2){
      p0[s,ss] ~ dunif(0,0.5)
      sigma[s,ss] ~ dunif(0,6)
    }#ss
  }#s
  

  for(i in 1:n.individuals){
      y[i,1:n.maxDets] ~ dbinomLocal_normal(
        size = size[1:n.detectors],
        p0 = p0[state[i],sex[i]+1],
        sigma = sigma[state[i],sex[i]+1],
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
##---- Set model constants, data & parameter simulated values (==inits)
simConstants <- list( n.individuals = 500,
                      n.habWindows = habitat$n.HabWindows,
                      n.states = 3,
                      n.detectors = detectors$n.detectors, 
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      n.maxDets = 40,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

simData <- list( lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords,
                 habIntensity = rep(1, habitat$n.HabWindows),
                 alpha = matrix(1,3,2),
                 size = detectors$trials,
                 detCoords = detectors$scaledCoords,
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects$habitatGrid)

simInits <- list( "psi" = 0.8,
                  "rho" = 0.5,
                  "theta" = matrix(c(0.5,0.45,0.05,
                                     0.5,0.3,0.2),3,2),
                  #"p0" = c(0.1,0.1,0.02),
                  #"sigma" = c(1,1,3))
                   "p0" = matrix(c(0.1,0.1,0.02,
                                   0.1,0.1,0.03),3,2),
                   "sigma" = matrix(c(1,1,2,
                                      1,1,3),3,2)) 

##---- Create the nimble model object
simModel <- nimbleModel( code = modelCode,
                         constants = simConstants,
                         inits = simInits,
                         data = simData,
                         check = FALSE,
                         calculate = FALSE) 

CsimModel <- compileNimble(simModel)

##---- Get nodes to simulate
nodesToSim <- CsimModel$getDependencies(c(names(simInits), "habIntensity"),
                                       self = F,
                                       includeData = F,
                                       downstream = T,
                                       returnScalarComponents = T)

nodesToSim.simple <- unique(sapply(strsplit(nodesToSim, "\\["), "[", 1))
nodesToSim.simple

##---- Simulate from NIMBLE 
CsimModel$simulate(nodesToSim)
CsimModel$calculate()
plot(habitat$loScaledCoords)
points( CsimModel$sxy[,1],
        CsimModel$sxy[,2],
        pch = 3, col = "red")

# CsimModel$simulate(nodesToSim)
# CsimModel$calculate()
# plot(habitat$loScaledCoords)
# points( CsimModel$sxy[,1],
#         CsimModel$sxy[,2],
#         pch = 3, col = "red")




## ------   3. BUNDLE NIMBLE DATA/INITIAL VALUES/... ------
nimConstants <- simConstants 
nimData <- simData
nimInits <- simInits

##---- Identify which individuals were detected
whichDet <- CsimModel$y[ ,1] > 0
whichNotDet <- CsimModel$y[ ,1] == 0

##---- Add simulated data 
nimData$y <- CsimModel$y

z <- z.init <- CsimModel$z
z[whichNotDet] <- NA
z.init[whichDet] <- NA
nimData$z <- z
nimInits$z <- z.init

sex <- sex.init <- CsimModel$sex
sex[whichNotDet] <- NA
sex.init[whichDet] <- NA
nimData$sex <- sex
nimInits$sex <- sex.init

state <- state.init <- CsimModel$state
state[whichNotDet] <- NA
state.init[whichDet] <- NA
nimData$state <- state
nimInits$state <- state.init

nimInits$sxy <- CsimModel$sxy

nimParams <- c("N", "p0", "sigma", "psi", "rho", "theta")
#"z", "sxy")




## ------   4. FIT MODEL -----
##---- Set-up NIMBLE model
nimModel <- nimbleModel( code = modelCode,
                         constants = nimConstants,
                         data = nimData,
                         inits = nimInits,
                         check = F,       
                         calculate = F)  
nimModel$calculate()

conf <- configureMCMC( model = nimModel,
                       monitors = nimParams,
                       thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
for(c in 1:3){
  runMCMCbites( mcmc = Cmcmc,
                bite.size = 1000,
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
