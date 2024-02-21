## -------------------------------------------------------------------------- ##
## ------------------------ ALPINE WOLF SC ---------------------------------- ##
## -------------------------------------------------------------------------- ##
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
library(tidyr)
library(lubridate)
library(stars)
library(units)
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(gridExtra)
library(MetBrewer)
library(data.table)
library(ggplot2)



## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
modelName = "AlpineWolf.SC.1.1.1"
thisDir <- file.path(analysisDir, modelName)

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}

## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 30000,
                dOK = 100, 
                dmax = 25)


## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
##---- Polygon of Italy and neighbouring countries
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

##--- Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() 

##---- create raster with desired resolution
grid.r <- st_as_stars(st_bbox(studyArea),
                      dx = habitat$resolution,
                      dy = habitat$resolution)
##---- mask and crop
grid.r <- grid.r[studyArea, ] 
##---- convert from stars to sf objects
grid <- st_as_sf(grid.r)
grid <- grid[studyArea, ] 
##---- give an ID to each grid cell 
grid$id <- 1:nrow(grid)
##---- Remove the "values" column
grid <- grid[-1]

##---- create centroids 
grid$centroids <- st_centroid(grid) %>% 
  st_geometry() 




## ------   2. CAMERA TRAPS ------
##---- Load GPS positions of Camera Traps
ct <- read_sf(file.path(dataDir,"GISData/CameraTraps/Ctraps/ct_alpi_240213.shp"))
sum(duplicated(ct[,c("coord_x","coord_y")]))

##---- Convert dates
ct$date_st <- parse_date_time(ct$date_st, orders = c('dmy'))
ct$date_en <- parse_date_time(ct$date_en, orders = c('dmy'))
ct$tot_attivi <- as.numeric(ct$tot_attivi) 

##---- Aggregate camera-traps with the same location 
# ct <- ct %>% group_by(coord_x,coord_y) %>%
#   summarize(date_st = min(date_st,na.rm = T),
#             date_en = max(date_en,na.rm = T),
#             non_workin = sum(`non-workin`, na.rm = T),
#             tot_attivi = sum(tot_attivi, na.rm = T) 
#   )

##---- Set-up new identifier for each camera-trap
ct$id <- paste0("FT", 1:nrow(ct))

##---- Plot check
plot(studyArea, col = "steelblue")
plot(ct$geometry, col = "blue", pch = 3, add = T)



## ------   3. PICTURES/SIGHTNG DATA ------
##---- Load all images available
pics <- read_sf(file.path(dataDir,"GISData/CameraTraps/photos/ft_alpi_photos_240213.shp"))
dim(pics)

##---- Use lubridate to clean dates
# The orders argument is a character vector containing the possible date-time parsing 
# formats in the order they should be tested. So by giving c('dmy', 'ymd'), lubridate 
# will try to parse all strings as Day, Month, Year format. If it can't do that 
# successfully (for example, the date 2021-03-18 won't work as there is no 2021th day), 
# it will try the next in the list until all strings are parsed, or all orders exhausted.
pics$date <- parse_date_time(pics$date, orders = c('dmy','mdy','ymd'))
pics$Year <- as.numeric(format(pics$date,"%Y"))
pics$Month <- as.numeric(format(pics$date,"%m"))
pics$n_wolves <- as.numeric(pics$n_wolves)

##---- Filter out samples outside the monitoring period
##---- Keep NAs (mostly from Val d'Aosta with weird format)
pics <- filter(pics, Month %in% c(1,2,3,4,10,11,12) | is.na(Month)) 
dim(pics)

##-- Filter out pictures of poor quality
pics <- filter(pics, C1 %in% 'C1')  
dim(pics)

##---- Create unique picture ID
pics$uniqueID <- 1:nrow(pics)

##---- Plot check
plot(pics$geometry, col = "black", add = T)



## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. DETECTORS ------
## ------     1.1. DETECTORS CHARACTERISTICS ------
detectors <- ct[ ,c("id","tot_attivi")]
detectors$coords <- st_coordinates(detectors)
dimnames(detectors$coords) <- list(1:nrow(detectors$coords),
                                   c("x","y"))
n.detectors <- nrow(detectors)



## ------   2. HABITAT ------
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (detector grid + buffer)
habitat$polygon <- st_buffer(st_union(grid),
                             habitat$buffer)

##---- Remove un-suitable habitat (e.g. seas)
habitat$polygon <- st_intersection( st_union(countries),
                                    habitat$polygon,
                                    drop_lower_td = TRUE)

##---- Create habitat raster (keep raster cells with > 50% habitat)
##---- (the huge buffer is used to make sure all models have comparable rasters)
habitat$raster <- raster(extent(st_bbox(st_buffer(st_union(studyAreaGrid),
                                                  100000))))

##---- Trick to get % of habitat
res(habitat$raster) <- habitat$resolution/10
habitat$raster <- fasterize( countries,
                             habitat$raster, )
habitat$raster[is.na(habitat$raster)]<- 0
habitat$raster <- aggregate(habitat$raster, fact = 10, sum, na.rm = TRUE)
habitat$raster[habitat$raster[ ] < 50] <- NA
habitat$raster[habitat$raster[ ] >= 50] <- 1
# plot(habitat$raster)


##---- Mask and crop habitat to the buffer
habitat$raster <- mask(habitat$raster, st_as_sf(habitat$polygon))
habitat$raster <- crop(habitat$raster, st_as_sf(habitat$polygon))

##---- Create habitat grid
habitat$grid <- st_as_sf(rasterToPolygons( habitat$raster,
                                           fun = function(x)x == 1))
habitat$grid$id <- 1:nrow(habitat$grid)
habitat$grid <- habitat$grid[,c(2,3)]
st_crs(habitat$grid) <- st_crs(habitat$polygon)

##---- Create Habitat matrix of cell ID
habitat$matrix <- habitat$raster
habitat$matrix[] <- 0

##----  Identify suitable habitat cells
isHab <- which(habitat$raster[]==1)

##---- Cell ID starts from the top left corner, increments left to right and
##---- up to down
habitat$matrix[isHab] <- 1:length(isHab)

##---- Convert to matrix
habitat$matrix <- as.matrix(habitat$matrix)
habitat$binary <- habitat$matrix  ##(required by the getLocalObject function; the function could be changed to deal with matrix others than 0/1)
habitat$binary[habitat$binary > 0] <- 1

##---- Obtain xy coordinates of habitat cells
habitat$coords <- coordinates(habitat$raster)[isHab, ]
dimnames(habitat$coords) <- list(1:length(isHab), c("x","y"))
habitat$sp <- SpatialPoints(coords = habitat$coords,
                            proj4string = crs(habitat$polygon))

##---- Retrieve habitat windows corners
habitat$lowerCoords <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)


##---- Visual plotting to check if everything is right
plot(habitat$raster)
plot(st_geometry(countries), add = T)
plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))
plot(ct,add=T, col="red")



## ------     2.2. HABITAT COVARIATES ------
## ------       2.2.1. HUMAN POPULATION DENSITY -------
##---- Load human population density data
POP_raw  <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
plot(POP_raw)

## Aggregate to the detector resolution
POP <- raster::aggregate( x = POP_raw ,
                          fact = habitat$resolution/res(POP_raw ),
                          fun = mean)
POP <- raster::focal(x = POP,
                     w = matrix(1,3,3),
                     fun = mean,
                     na.rm = T)
plot(POP)
plot(st_geometry(grid), add = T)

## Extract scaled human pop density
habitat$grid$pop <- POP %>%
  raster::extract(y = habitat$sp) %>%
  scale()

## Extract log(human pop density + 1)
log_POP <- POP
log_POP[ ] <- log(POP[]+1)

par(mfrow=c(1,2))
plot(POP)
plot(st_geometry(countries),add=T)
plot(log_POP)
plot(st_geometry(countries),add=T)

habitat$grid$log_pop <- log_POP %>%
  raster::extract(y = habitat$sp) %>%
  scale()


## ------       2.2.2. CORINE LAND COVER ------
##---- Load Corine Land Cover data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))

## Legend :
## 1	= developed; 2 = agriculture; 3 = forest; 4 = herbaceous; 5 =	bare rock;
## 6 = sea features; 7 = inland water; 9 = perpetual snow
layer.names <- c("developed","agriculture","forest","herbaceous","bare rock",
                 "sea features","inland water","NA","perpetual snow")
layer.index <- c(1,2,3,4,5,9)
par(mfrow = c(1,2))

CLC.df <- as.data.frame(matrix(NA, habitat$n.HabWindows, length(layer.index)))
names(CLC.df) <- layer.names[layer.index]
CLC.df$id <- 1:habitat$n.HabWindows

COV <- list()
for(l in 1:length(layer.index)){
  temp <- CLC
  temp[!temp[] %in% layer.index[l]] <- 0
  temp[temp[] %in% layer.index[l]] <- 1
  plot(temp)
  plot(st_geometry(habitat$grid), add = T)
  
  COV[[l]] <- raster::aggregate( x = temp,
                                 fact = habitat$resolution/res(temp),
                                 fun = mean)
  
  COV[[l]] <- raster::focal(x = COV[[l]],
                            w = matrix(1,3,3),
                            fun = mean,
                            na.rm = T)
  
  plot(COV[[l]])
  plot(st_geometry(habitat$grid), add = T)
  CLC.df[ ,l] <- raster::extract( x = COV[[l]],
                                  y = habitat$sp)
  print(layer.index[l])
}#c

## Store in habitat grid
habitat$grid <- left_join(habitat$grid, CLC.df, by = "id")



## ------       2.2.3. IUCN PRESENCE ------
##---- Load IUCN presence data
list.files(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid"))
iucn_2012_1 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/Clip_2012_12_01_Wolves_permanent.shp"))
iucn_2012_1$SPOIS <- 3
iucn_2012_1 <- st_transform(iucn_2012_1, st_crs(studyArea))
plot(st_geometry(studyArea))
plot(st_geometry(iucn_2012_1), add = T, col = "lightskyblue4")

iucn_2012_2 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/Clip_2012_12_01_Wolves_sporadic.shp"))
iucn_2012_2$SPOIS <- 1
iucn_2012_2 <- st_transform(iucn_2012_2, st_crs(studyArea))
plot(st_geometry(iucn_2012_2), add = T, col = "lightskyblue2")

iucn_2018 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/2018_06_06_Wolf_IUCN_Redlist.shp"))
iucn_2018 <- st_transform(iucn_2018, st_crs(studyArea))
plot(st_geometry(studyArea))
plot(iucn_2018[ ,"SPOIS"], add = T)
## Code back to numeric
iucn_2018$SPOIS <- ifelse(iucn_2018$SPOIS == "Sporadic", 1, 3)

## Extract LCIE wolf permanent presence in each habitat grid cell
intersection <- st_intersection(habitat$grid, iucn_2012_1) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(IUCN = sum(iucn*SPOIS)/2.5e+07)

habitat$grid <- habitat$grid %>%
  left_join(intersection, by = "id")
habitat$grid$IUCN[is.na(habitat$grid$IUCN)] <- 0
plot(habitat$grid[,"IUCN"])

## Extract LCIE wolf sporadic presence in each habitat grid cell
intersection <- st_intersection(habitat$grid, iucn_2012_2) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/2.5e+07)

tmp <- habitat$grid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitat$grid$IUCN <- habitat$grid$IUCN + tmp$iucn_2
plot(habitat$grid[,"IUCN"])

## Extract LCIE wolf presence in each habitat grid cell
intersection <- st_intersection(habitat$grid, iucn_2018) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/2.5e+07)

tmp <- habitat$grid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitat$grid$IUCN <- habitat$grid$IUCN + tmp$iucn_2
plot(habitat$grid[,"IUCN"])

habitat$grid$IUCN <- scale(habitat$grid$IUCN)



## ------   3. RESCALE HABITAT & DETECTORS ------
##---- Rescale habitat and detector coordinates
scaledCoords <- scaleCoordsToHabitatGrid(
  coordsData = detectors$coords,
  coordsHabitatGridCenter = habitat$coords)

detectors$scaledCoords <- scaledCoords$coordsDataScaled
habitat$scaledCoords <- scaledCoords$coordsHabitatGridCenterScaled
habitat$loScaledCoords <- habitat$scaledCoords - 0.5
habitat$upScaledCoords <- habitat$scaledCoords + 0.5



## ------   4. IDENTIFY LOCAL CAMERA-TRAPS ------
localObjects <- getLocalObjects(
  habitatMask = habitat$binary,
  coords = scaledCoords$coordsDataScaled,
  dmax = habitat$dmax)



## ------   5. DETECTION DATA ------
## ------     5.1. ASSIGN PICTURES TO CAMERA-TRAPS ------
##---- Calculate distance between detections and sub-detectors
closest <- nn2( st_coordinates(detectors),
                st_coordinates(pics),
                k = 1,
                searchtype = "radius",
                radius = 30000)


##---- Identify and plot pictures to check
table(closest$nn.dists)
toCheck <- pics[closest$nn.dists > 0, ]

plot(studyArea, col = "steelblue")
plot(ct$geometry, col = "blue", pch = 3, add = T)
plot(pics$geometry, col = "black", add = T)
plot(pics$geometry[toCheck], col = "red",pch=19, cex = 0.8, add = T)


##---- Assign each picture to a camera-trap based on minimum distance
pics$CT_id <- detectors$id[closest$nn.idx] 

##---- Remove pictures further than 'dOK' from a registered camera trap
pics <- pics[closest$nn.dists <= habitat$dOK, ]




## ------     5.2. AGGREGATE DATA (NOT FOR NOW) ------
# ##---- Aggregate pictures per DAY 
# pics_d <- pics %>% 
#   group_by(day = lubridate::day(date), id) %>% 
#   summarise(n_wolves = max(n_wolves))
# hist(pics_d$n_wolves)
# dim(pics_d)
# 
# ##---- Aggregate pictures per WEEK 
# pics_w <- pics %>% 
#   group_by(week = lubridate::week(date), id) %>% 
#   summarise(n_wolves = max(n_wolves))
# hist(pics_w$n_wolves)
# dim(pics_w)
# 
# ##---- Aggregate pictures per MONTH 
# pics_m <- pics %>% 
#   group_by(Month, id) %>% 
#   summarise(n_wolves = max(n_wolves))
# hist(pics_m$n_wolves)
# dim(pics_m)




## ------     5.3. DETECTION VECTOR : y ------
##-- Create a vector of camera-trap specific observations
y <- y.visits <- rep(0, n.detectors)

##-- Loop over camera-traps
for(j in 1:nrow(detectors)){
  ##-- Identify this detector
  thisCT <- detectors$id[j]
  ##-- Subset pictures to this camera-trap
  thesePics <- filter(pics, CT_id == thisCT)
  if(nrow(thesePics)>0){
    ##-- Calculate total number of wolves detected at this camera-trap
    y[j] <- sum(thesePics$n_wolves)
    ##-- Calculate total number of pictures (=events) at this camera-trap
    y.visits[j] <- nrow(thesePics)
  }#if
}#j





## -----------------------------------------------------------------------------
## ------ III. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] ~ dnorm(0.0,0.01)
  }#c
  
  ##-- Intensity of the AC distribution point process
  habIntensity[1:n.habWindows] <- exp( hab.covs[1:n.habWindows,1:n.habCovs] %*% betaHab[1:n.habCovs])
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
  sigma ~ dunif(0,10)
  lambda0 ~ dunif(0,10)

  for(i in 1:M){ 
    lambda[i,1:J] <- calculateLocalLambda(
      lambda0 = lambda0,
      sigma = sigma,
      s = s[i,1:2],
      trapCoords = trapCoords[1:J,1:2],
      localTrapsIndices = localTrapsIndices[1:n.habWindows,1:n.localIndicesMax],
      localTrapsNum = localTrapsNum[1:n.habWindows],
      resizeFactor = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      indicator = z[i])
  }#i
  
  for(j in 1:J){
    bigLambda[j] <- sum(lambda[1:M,j]) 
    y[j] ~ dpois(bigLambda[j]*oper[j])
  }#j
  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
  D <- N/area
})



## ------   2. BUNDLE DATA ------
##---- Set model data
area <- st_area(grid) %>%
  drop_units() %>%
  sum()

M <- 3000

nimData <- list( y = y,
                 area = area,
                 oper = detectors$tot_attivi,
                 habitatGrid = localObjects$habitatGrid,
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 trapCoords = detectors$scaledCoords,
                 lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords,
                 hab.covs = cbind.data.frame(
                   "bare rock" = habitat$grid$`bare rock`,
                   "herbaceous" = habitat$grid$`herbaceous`,
                   "forest" = habitat$grid$`forest`,
                   "pop" = habitat$grid$pop,
                   "IUCN" = habitat$grid$`IUCN`))

##---- Set model constants
nimConstants <- list( M = M,
                      J = length(y),
                      n.habWindows = habitat$n.HabWindows,
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2],
                      n.habCovs = ncol(nimData$hab.covs))

##---- Set parameters to track 
nimParams <- c("N", "D", "lambda0", "sigma", "psi", "betaHab")
nimParams2 <- c("z", "s")




## ------   3. NIMBLE INITIAL VALUES ------
##---- Generate initial individual AC locations
sInits <- matrix(NA, nrow = M, ncol = 2)
for (i in 1:M) {
  sInits[i, ] <- rbernppAC( n = 1,
                            lowerCoords = nimData$lowerHabCoords,
                            upperCoords = nimData$upperHabCoords,
                            logIntensities = log(rep(1,habitat$n.HabWindows)),
                            logSumIntensity = log(habitat$n.HabWindows),
                            habitatGrid = nimData$habitatGrid,
                            numGridRows = nimConstants$y.max,
                            numGridCols = nimConstants$x.max)
}#i

##---- Generate initial number of days operated (only when NA)
operInits <- rpois( n = nrow(detectors), mean(nimData$oper, na.rm = T))
operInits[!is.na(nimData$oper)] <- NA

##---- Create a list of random initial values (one set per chain)
nimInits.list <- list()
for(c in 1:4){
  nimInits.list[[c]] <- list( 
    sigma = 0.7,
    lambda0 = runif(1,0,2),
    oper = operInits,
    s = sInits,
    lambda.oper = 1,
    z = rbinom(M,1,0.6),
    psi = 0.6,
    betaHab = rep(0,nimConstants$n.habCovs),
    habIntensity = rep(1,habitat$n.HabWindows))
}#c




## ------   4. SAVE NIMBLE INPUT -----
for(c in 1:4){
  nimInits <- nimInits.list[[c]]
  save( modelCode,
        nimData,
        nimConstants,
        nimInits,
        nimParams,
        nimParams2,
        file = file.path(thisDir, "input",
                         paste0(modelName, "_", c, ".RData")))
}#c




## -----------------------------------------------------------------------------
## ------ IV. FIT NIMBLE MODEL -----
for(c in 1:4){
  load( file.path(thisDir, "input", paste0(modelName,"_", c ,".RData")))
  
  ##---- Create the nimble model object
  nimModel <- nimbleModel( code = modelCode,
                           constants = nimConstants,
                           inits = nimInits,
                           data = nimData,
                           check = FALSE,
                           calculate = FALSE) 
  nimModel$calculate()
  
  
  
  
  ##---- Compile the nimble model object to C++
  CsimModel <- compileNimble(nimModel)
  CsimModel$calculate()
  
  ##---- Configure and compile the MCMC object 
  conf <- configureMCMC( model = nimModel,
                         monitors = nimParams,
                         thin = 10,
                         monitors2 = nimParams2,
                         thin2 = 40)
  
  
  Rmcmc <- buildMCMC(conf)
  compiledList <- compileNimble(list(model = nimModel,
                                     mcmc = Rmcmc))
  Cmodel <- compiledList$model
  Cmcmc <- compiledList$mcmc
  
  ##---- Run nimble MCMC in multiple bites
  mcmcRuntime <- system.time(
    runMCMCbites( mcmc = Cmcmc,
                  model = Cmodel,
                  bite.size = 100,
                  bite.number = 1,
                  path = file.path(thisDir, "output", modelName, "_", c),
                  save.rds = TRUE))  
  
  print(mcmcRuntime)
}#c 



## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 0)

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput)
graphics.off()


##---- Process and save MCMC samples
res <- ProcessCodaOutput(nimOutput$samples)
res_sxy <- ProcessCodaOutput(nimOutput$samples2)
save(res, res_sxy, 
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))




## ------   1. CALCULATE AC-BASED DENSITY ------
##---- Load processed MCMC samples
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
load(file.path(thisDir, "input", paste0(modelName, "_1.RData")))

##----  Extract the habitat raster
habitat.r <- habitat$raster

##---- Set cells outside the habitat to NA
habitat.r[habitat.r == 0] <- NA

##---- Create a matrix of raster cellsIDs
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = dim(habitat.r)[1],
                      ncol = dim(habitat.r)[2],
                      byrow = TRUE)

##---- Create a matrix of admin units binary indicators
##---- (rows == regions ; columns == habitat raster cells)
hab.rgmx <- rbind(habitat$raster[] == 1)
hab.rgmx[is.na(hab.rgmx)] <- 0
row.names(hab.rgmx) <- "habitat"

##---- Calculate density
WA_Density <- GetDensity(
  sx = res_sxy$sims.list$s[ , ,1],
  sy = res_sxy$sims.list$s[ , ,2],
  z = res_sxy$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = hab.rgmx)


##---- Create a matrix of italy 
##---- (rows == regions ; columns == habitat raster cells)
regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = 0)
regions.r <- regions.r + habitat$Italia - 1
plot(regions.r)
table(regions.r[], useNA = "always")

regions.unique <- na.omit(unique(regions.r[]))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions$DEN_UTS[regions.unique] 

alps.r <- fasterize(sf = st_as_sf(alps),
                    raster = habitat.r,
                    background = 0)
alps.r <- alps.r + habitat$Italia - 1
plot(alps.r)
table(alps.r[], useNA = "always")
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"


##---- Calculate overall density
WA_Italy <- GetDensity(
  sx = res_sxy$sims.list$s[ , ,1],
  sy = res_sxy$sims.list$s[ , ,2],
  z = res_sxy$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = regions.rgmx)

WA_status <- list()
for(s in 0:1){
  WA_status[[s+1]] <- list()
  for(ss in 1:3){
    status <- (res_sxy$sims.list$z == 1) &
      (res_sxy$sims.list$sex == s) &
      (res_sxy$sims.list$status == ss)
    
    WA_status[[s+1]][[ss]] <- GetDensity(
      sx = res_sxy$sims.list$s[ , ,1],
      sy = res_sxy$sims.list$s[ , ,2],
      z = status,
      IDmx = habitat.id,
      aliveStates = 1,
      returnPosteriorCells = T,
      regionID = regions.rgmx)
  }
}

##---- Create a matrix of admin units binary indicators
##---- (rows == regions ; columns == habitat raster cells)
comp.rgmx <- rbind(habitat$comparison[] == 1)
comp.rgmx[is.na(comp.rgmx)] <- 0
row.names(comp.rgmx) <- "comparison"

##---- Calculate density
WA_Comp <- GetDensity(
  sx = res_sxy$sims.list$s[ , ,1],
  sy = res_sxy$sims.list$s[ , ,2],
  z = res_sxy$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = comp.rgmx)


##---- Create a matrix of Wolf Presence grid 
##---- (rows == regions ; columns == habitat raster cells)
regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = 0)
regions.r <- regions.r + habitat$extraction - 1
plot(regions.r)
table(regions.r[], useNA = "always")

regions.unique <- na.omit(unique(regions.r[]))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions$DEN_UTS[regions.unique] 

alps.r <- fasterize(sf = st_as_sf(alps),
                    raster = habitat.r,
                    background = 0)
alps.r <- alps.r + habitat$extraction - 1
plot(alps.r)
table(alps.r[], useNA = "always")
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"

##---- Calculate density
WA_Extract <- GetDensity(
  sx = res_sxy$sims.list$s[ , ,1],
  sy = res_sxy$sims.list$s[ , ,2],
  z = res_sxy$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = regions.rgmx)




## ------   2. DENSITY ------
pdf(file = file.path(thisDir, paste0(modelName,"_Density.pdf")),
    width = 20, height = 15)

## ------     2.1. DENSITY MAP ------
##---- Set color scale
maxDens <- max(WA_Density$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)
par(mfrow = c(2,2), mar = c(4,4,0,0))

##---- Plot overall density raster
meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_Density$MeanCell
meanDensity.R[is.na(habitat.r[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot( st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Density$summary["Total",1],1),
                     " [", round(WA_Density$summary["Total",4],1), " ; ",
                     round(WA_Density$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot Italian density raster
ital.R <- habitat$Italia
ital.R[] <- WA_Italy$MeanCell
ital.R[is.na(habitat$Italia[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( ital.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot( st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Italy$summary["Total",1],1),
                     " [", round(WA_Italy$summary["Total",4],1), " ; ",
                     round(WA_Italy$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot density raster for comparison between models
comp.R <- habitat.r
comp.R[ ] <- WA_Comp$MeanCell
comp.R[is.na(habitat$comparison[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( comp.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot(st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Comp$summary["Total",1],1),
                     " [", round(WA_Comp$summary["Total",4],1), " ; ",
                     round(WA_Comp$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot density raster for extraction
extract.R <- habitat.r
extract.R[ ] <- WA_Extract$MeanCell
extract.R[is.na(habitat$extraction[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( extract.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot(st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Extract$summary["Total",1],1),
                     " [", round(WA_Extract$summary["Total",4],1), " ; ",
                     round(WA_Extract$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


## ------  REALIZED vs. PREDICTED DENSITY  
##-- Calculate relative density
relativeDens.r <- meanDensity.R
relativeDens.r[] <- meanDensity.R[]/sum(meanDensity.R[], na.rm = T)

##-- Calculate relative point-process intensity
intensity.r <- habitat.r
intensity.r[intensity.r[] > 0] <- c(exp(res$mean$betaHab %*% t(nimData$hab.covs)))
relativeInt.r <- intensity.r
relativeInt.r[] <- relativeInt.r[]/sum(relativeInt.r[],na.rm=T)

##-- Set color scale
maxDens <- max(c(relativeInt.r[], relativeDens.r[]), na.rm = T)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)

par(mfrow = c(1,2), mar = c(1,1,1,1))
plot( habitat$polygon, border = grey(0.3))
mtext( text = "Realized density", side = 3, line = -4, font = 2)
plot( relativeDens.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T)
plot( st_geometry(countries), add = T, lwd = 2)

plot( habitat$polygon, border = grey(0.3))
mtext( text = "Predicted density", side = 3, line = -4, font = 2)
plot( relativeInt.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T)
plot( st_geometry(countries), add = T, lwd = 2)

## Add density legend
plot( relativeInt.r, breaks = cuts, col = col,
      legend.width = 1, legend.only = TRUE,
      axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 5),
                        labels = round(seq(0,maxDens,length.out = 5)*100, digits = 2), 
                        cex.axis = 0.6),
      legend.args = list( text = 'Relative density (% of pop.km2)', line = -2,
                          side = 4, font = 2, cex = 0.8))



## ------  DETECTION CENTROIDS
par(mfrow = c(1,1), mar = c(4,4,4,4))

##-- Set color scale
maxDens <- max(meanDensity.R[], na.rm = T)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)

centroids <- st_drop_geometry(ngs) %>%
  group_by(Genotype.ID) %>%
  summarise(X = mean(CoordX),               ## Get total length searched in each detector grid cell
            Y =  mean(CoordY)) 
coordinates(centroids) <- cbind.data.frame(centroids$X,
                                           centroids$Y)
proj4string(centroids) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
centroids <- st_as_sf(centroids)

plot( habitat$polygon, border = grey(0.3))
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot(  habitat$polygon, add = T)
plot(st_geometry(countries), add = T, lwd=2)
plot(centroids,
     add = T, pch = 21, cex = 1.2,
     bg = adjustcolor("red",0.2),
     col = "red")
mtext( text = "Wolf density with\ndetected individuals' centroid",
       side = 3, font = 2, cex = 2)


## ------  REGIONAL ABUNDANCES
##-- Export as .pdf
par(mfrow = c(1,1))
plot.new()
grid.table(round(WA_Italy$summary,1))
mtext(text = "Abundance estimates\nby region",
      side = 3, outer = T, line = -20, cex = 2, font = 2)

sex <- c("female","male")
status <- c("alpha","pups","others")
par(mfrow = c(2,3))
for(s in 1:length(sex)){
  for(ss in 1:length(status)){
    par(mfrow = c(1,1))
    plot.new()
    grid.table(round(WA_status[[s]][[ss]]$summary,1))
    mtext(text = paste0("Abundance estimates\nfor ", sex[s], "_", status[ss]),
          side = 3, outer = T, line = -10, font = 2)
  }
}



## ------     2.2. DENSITY EFFECT PLOT ------
# par(mfrow = c(2,2), mar = c(6,6,0,0))
# covNames <- names(nimData$hab.covs)# c("alpine", "forest", "IUCN")
# pred.hab.covs <- apply(nimData$hab.covs,
#                        2,
#                        function(c){
#                          seq(min(c),
#                              max(c),
#                              length.out = 100)})
# mean.hab.covs <- apply(nimData$hab.covs, 2, mean)
# cols <- hcl.colors(length(mean.hab.covs))
# 
# for(b in 1:ncol(res$sims.list$betaHab)){
#   intensity <- do.call( rbind,
#                         lapply(1:length(res$sims.list$betaHab[ ,b]),
#                                function(x){
#                                  exp(res$sims.list$betaHab[x,b] * pred.hab.covs[ ,b] + 
#                                        sum(res$sims.list$betaHab[x,-b] * mean.hab.covs[-b]))
#                                }))
#   
#   
#   mean.int <- colMeans(intensity)
#   quant.int <- apply(intensity, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
#   maxD <- round(max(quant.int),3)
#   plot( x = pred.hab.covs[ ,b],
#         y = quant.int[2, ],
#         type = "n", ylim = c(0, maxD), xlim = range(pred.hab.covs[ ,b]),
#         ylab = "Density", xlab = covNames[b], axes = FALSE)
#   minCov <- min(st_drop_geometry(habitat$grid[ ,covNames[b]]))
#   maxCov <- max(st_drop_geometry(habitat$grid[ ,covNames[b]]))
#   xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
#   axis(1,
#        at = round(seq(min(pred.hab.covs[ ,b]), max(pred.hab.covs[ ,b]), length.out = 10),3),
#        labels = xLabels, cex = 2,
#        tck = 0.01, las = 1, hadj = 0.5)
#   axis(2, at = seq(0,maxD,length.out = 6),
#        labels = seq(0,maxD,length.out = 6),
#        tck = 0.01, las = 1, cex = 2)
#   
#   polygon(x = c(pred.hab.covs[ ,b],rev(pred.hab.covs[ ,b])),
#           y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
#           col = adjustcolor(cols[b], alpha.f = 0.5))
#   
#   points( x = pred.hab.covs[ ,b],
#           y = quant.int[2,],
#           lwd = 2, type = "l", col = cols[b])
#   
# }
# graphics.off()
# 
# 
# 
# 
# 
## ------     2.2. RJ-MCMC PLOTS ------

##---- Get model features
n.chains <- length(nimOutput$samples)
n.iterations <- dim(nimOutput$samples[[1]])[1]

##---- Get target covariates
covNames <- colnames(nimData$hab.covs)

zRJ.wide <- data.table(res$sims.list$zRJ)
dimnames(zRJ.wide) <- list(NULL, covNames)

betas.wide <- data.table(res$sims.list$betaHab)
dimnames(betas.wide) <- list(NULL, covNames)


##---- List model combinations
mods <- apply(zRJ.wide, 1, function(x){paste(covNames[x == 1], collapse = "+")})

betas.wide$model <- zRJ.wide$model <- gsub("\\(Intercept\\)\\+", "", mods)
betas.wide$chain <- zRJ.wide$chain <- rep(1:n.chains, each = n.iterations)
betas.wide$iteration <- zRJ.wide$iteration <- rep(1:n.iterations, n.chains)

zRJ.df <- melt(zRJ.wide, id.vars = c("iteration", "chain", "model"))
names(zRJ.df) <- c("iteration", "chain", "model", "variable", "value")

betas.df <-  melt(betas.wide, id.vars = c("iteration", "chain", "model"))
names(betas.df) <-  c("iteration", "chain", "model", "variable", "value")

betas.df$value[zRJ.df$value == 0] <- NA

betas.aggr <- betas.df %>%
  group_by(variable) %>%
  summarise(p.inclusion = mean(!is.na(value)))

betas.df <- merge(betas.df, betas.aggr)

included <- zRJ.df$value == 1

betas.df <- betas.df[included, ]
zRJ.df <- zRJ.df[included, ]

betas.df <- betas.df[order(betas.df$variable, betas.df$model, betas.df$chain), ]
zRJ.df <- zRJ.df[order(zRJ.df$variable, zRJ.df$model, zRJ.df$chain), ]

myfun1 <- function(x) 1:length(x)

temp <- betas.df %>% group_by(variable, model, chain) %>%
  summarize(iteration.model = myfun1(value))

betas.df$iteration.model <- temp$iteration.model

aggr <- data.frame(table(betas.df$model) / length(betas.df$model))
names(aggr) <- c("model", "weight")
betas.df <- merge(betas.df, aggr)


##---- MODEL TALLY
aggr <- aggr[order(aggr$weight, decreasing = TRUE),]
aggr$model <- factor(aggr$model, levels = aggr$model)
ggplot(data = aggr,
       mapping =  aes(x = model, y = weight, alpha = weight)) +
  geom_col(fill =
             "magenta") + theme(axis.text.x = element_text(
               angle = 45,
               vjust = 1,
               hjust = 1
             )) + ylab("Weight") + xlab("Models")


##---- COEFFICIENT TRACE PLOTS (OVERALL)
ggplot(data = betas.df, aes(
  x = iteration,
  y = value,
  color = factor(chain))) +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  xlab("Iteration") +  theme(legend.position = "none")


##---- COEFFICIENT TRACE PLOTS (MODEL-SPECIFIC)
ggplot(data = betas.df, aes(
  x = iteration.model,
  y = value,
  color = factor(chain))) +
  geom_line() +
  facet_grid(variable ~ model, margins = FALSE, scales = "free") +
  xlab("Iteration") +  theme(legend.position = "none")


##---- PLOT COEFFICIENT ESTIMATES (OVERALL)
ggplot(betas.df, aes(value, variable, alpha = p.inclusion)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "turquoise",
    color = grey(1)) +
  geom_vline(xintercept = 0)


##---- PLOT COEFFICIENT ESTIMATES (MODEL-SPECIFIC)
ggplot(betas.df, aes(value, variable, alpha = weight)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "magenta",
    color = grey(1)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ model)





## ------     2.3. DENSITY EFFECT PLOT ------
par(mfrow = c(2,2), mar = c(6,6,0,0))
covNames <- names(nimData$hab.covs)# c("alpine", "forest", "IUCN")
pred.hab.covs <- apply(nimData$hab.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.hab.covs <- apply(nimData$hab.covs, 2, mean)
cols <- hcl.colors(length(mean.hab.covs))

for(b in 1:ncol(res$sims.list$betaHab)){
  ##-- Calculate intensity
  intensity <- do.call( rbind,
                        lapply(1:length(res$sims.list$betaHab[ ,b]),
                               function(x){
                                 exp(res$sims.list$betaHab[x,b] * pred.hab.covs[ ,b] + 
                                       sum(res$sims.list$betaHab[x,-b] * mean.hab.covs[-b]))
                               }))
  quant.int <- apply(intensity, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
  maxD <- round(max(quant.int),3)
  plot( x = pred.hab.covs[ ,b],
        y = quant.int[2, ],
        type = "n", ylim = c(0, maxD),
        xlim = range(pred.hab.covs[ ,b]),
        ylab = "Density",
        xlab = covNames[b],
        axes = FALSE)
  minCov <- min(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  maxCov <- max(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
  axis(1,
       at = round(seq(min(pred.hab.covs[ ,b]), max(pred.hab.covs[ ,b]), length.out = 10), 3),
       labels = xLabels, cex = 2,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxD,length.out = 6),
       labels = seq(0,maxD,length.out = 6),
       tck = 0.01, las = 1, cex = 2)
  
  
  # for(m in 1:nrow(aggr)){
  m <- 1
  ##-- Identify iterations for this model
  tmp <- as.data.frame(betas.wide[betas.wide$model == aggr$model[m],1:5])
  
  ##-- Calculate intensity
  intensity <- do.call( rbind,
                        lapply(1:nrow(tmp),
                               function(x){
                                 exp(tmp[x,b] * pred.hab.covs[ ,b] + 
                                       sum(tmp[x,-b] * mean.hab.covs[-b]))
                               }))
  quant.int <- apply(intensity, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
  polygon(x = c(pred.hab.covs[ ,b],rev(pred.hab.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = round(aggr$weight[m],2)))
  
  points( x = pred.hab.covs[ ,b],
          y = quant.int[2,],
          lwd = 2, type = "l", col = cols[b])
  
  #}# m
}# b
graphics.off()









## ------   3. DETECTION ------
pdf(file = file.path(thisDir, paste0(modelName,"_Detection.pdf")),
    width = 15, height = 15)

## ------     3.1. DETECTION MAP ------
sex <- c("female","male")
status <- c("alpha","pup","other")
for(ss in 1:2){
  for(s in 1:3){
    detectors$grid[ ,paste0("p0_",sex[ss],"_",status[s])] <- c(ilogit(logit(res$mean$p0[s,ss]) +
                                                                        res$mean$betaDet %*% t(nimData$det.covs)))
    # plot(st_geometry(detectors$grid[ ,paste0("p0.",ss,".",s)]),
    #      main = paste0("p0_",sex[ss],"_",status[s]))
    # plot(detectors$grid[ ,paste0("p0.",ss,".",s)],
    #      add = T,
    #      breaks = seq(0,1,0.05))    
  }
}
plot(detectors$grid[ ,c("p0_female_alpha", 
                        "p0_male_alpha",
                        "p0_female_pup",
                        "p0_male_pup",
                        "p0_female_other",
                        "p0_male_other")],
     key.pos = 4, breaks = seq(0,1,0.05)) 


## ------     3.2. DETECTION EFFECT PLOT ------
par(mfrow = c(3,2), mar = c(8,8,0,4))
covNames <- c("transect_L",
              "transect_qi",
              "snow_fall",
              "zone",
              "log_pop")

pred.det.covs <- apply(nimData$det.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.det.covs <- apply(nimData$det.covs, 2, mean)
cols <- hcl.colors(length(mean.det.covs))

for(b in 1:ncol(res$sims.list$betaDet)){
  p0 <- do.call( rbind,
                 lapply(1:length(res$sims.list$betaDet[,b]),
                        function(x){
                          ilogit(logit(res$sims.list$p0[x,1,1]) +
                                   res$sims.list$betaDet[x,b] * pred.det.covs[ ,b]+ 
                                   sum(res$sims.list$betaDet[x,-b] * mean.det.covs[-b]))
                        }))
  quant.int <- apply(p0, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
  maxp0 <- round(max(quant.int),3)
  plot( x = pred.det.covs[ ,b],
        y = quant.int[2, ],
        type = "n", ylim = c(0, maxp0), xlim = range(pred.det.covs[ ,b]),
        ylab = "p0", xlab = covNames[b], axes = FALSE)
  minCov <- min(st_drop_geometry(detectors$grid[ ,covNames[b]]))
  maxCov <- max(st_drop_geometry(detectors$grid[ ,covNames[b]]))
  xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
  axis(1,
       at = round(seq(min(pred.det.covs[ ,b]), max(pred.det.covs[ ,b]), length.out = 10),3),
       labels = xLabels, cex = 3,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxp0,length.out = 6),
       labels = seq(0,maxp0,length.out = 6),
       tck = 0.01, las = 1, cex = 3)
  
  polygon(x = c(pred.det.covs[ ,b],rev(pred.det.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = 0.5))
  
  points( x = pred.det.covs[ ,b],
          y = quant.int[2,],
          lwd = 2, type = "l", col = cols[b])
  
}
graphics.off()




## ------   4. PARAMETER TABLES ------
## ------     4.1. SCR PARAMETERS ------
#paramSimple <- sapply(strsplit(colnames(res$sims.list), split = '\\['), '[', 1)
params.simple <- names(res$mean)[!names(res$mean) %in% c("s","z")]

params.means <- do.call(c, lapply(params.simple, function(x)res$mean[[x]]))
params.sd <- do.call(c, lapply(params.simple, function(x)res$sd[[x]]))
params.q2.5 <- do.call(c, lapply(params.simple, function(x)res$q2.5[[x]]))
params.q97.5 <- do.call(c, lapply(params.simple, function(x)res$q97.5[[x]]))
params.Rhat <- do.call(c, lapply(params.simple, function(x)res$Rhat[[x]]))
params.n.eff <- do.call(c, lapply(params.simple, function(x)res$n.eff[[x]]))

params.summary <- cbind.data.frame( params.means, params.sd,
                                    params.q2.5, params.q97.5, params.Rhat,
                                    params.n.eff)

names(params.summary) <- c("mean", "sd", "2.5%CI", "97.5%CI", "Rhat", "n.eff")

## WRITE TABLE
# write.csv( x = params.summary,
#            file = file.path(myVars$WD, myVars$modelName, "TABLES",
#                             paste( myVars$modelName, "_params.csv", sep = "")))
# 
# print( xtable(params.summary, type = "latex"),
#        floating = FALSE,# scalebox=.8,
#        add.to.row = list(list(seq(1, nrow(params.summary), by = 2)),"\\rowcolor[gray]{.95} "),
#        file = file.path(myVars$WD, myVars$modelName,"TABLES",
#                         paste( myVars$modelName, "_params.tex", sep = "")))




## ------     4.2. ABUNDANCES ------
WA_Italy$summary



## ------   5. SPACE-USE ------
# ## Save input files for a few iterations
# iter <- seq(1, dim(res$sims.list$s)[1], by = 500)
# 
# WA_SpaceUse <- GetSpaceUse(  sx = res$sims.list$s[iter, ,1],
#                              sy = res$sims.list$s[iter, ,2],
#                              z = res$sims.list$z[iter, ],
#                              sigma = matrix(res$sims.list$sigma[iter,1,1],
#                                             length(iter),
#                                             ncol(res$sims.list$z)),
#                              habitatxy = habitat$scaledCoords,
#                              aliveStates = 1,
#                              regionID = it.rgmx,
#                              returnPosteriorCells = F)
# WA_SpaceUse$summary
# 
# ##-- Set color scale
# maxDens <- max(c(WA_Density$MeanCell, WA_SpaceUse$MeanCell))
# cuts <- seq(0, maxDens, length.out = 100)
# colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
# col <- colFunc(100)
# 
#
# ## ------       5.1. AC-based DENSITY MAP 
# ital.R <- habitat.r
# ital.R[ ] <- WA_Density$MeanCell
# ital.R[is.na(habitat$Italia[])] <- NA
# 
# plot( habitat$polygon, col = "gray80", border = grey(0.3))
# plot( ital.R, add = T,
#       breaks = cuts, col = col,
#       axes = F, box = F, bty = "n", legend = F)
# plot( habitat$polygon, add = T, border = grey(0.3))
# plot( st_geometry(countries), add = T, lwd = 2)
# 
# mtext( text = paste( "N = ", round(WA_Italy$summary["Total",1],1),
#                      " [", round(WA_Italy$summary["Total",4],1), " ; ",
#                      round(WA_Italy$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# 
# plot( meanDensity.R, breaks = cuts, col = col,
#       legend.width = 1, legend.only = TRUE,
#       axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
#                         labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
#                         cex.axis = 0.6),
#       legend.args = list( text = 'Density (ids.km2)', line = -2,
#                           side = 4, font = 2, cex = 0.8))
# 
# 
# 
# 
# ## ------       5.2. UD-based DENSITY MAP 
# habitat.r <- Habitat.list$habitat.r
# habitat.r[habitat.r == 0] <- NA
# admin.r <- rasterize(Cleanpoly, habitat.r)
# 
# meanDensity.R <- habitat.r
# meanDensity.R[meanDensity.R[]==1] <- WA_SpaceUse$MeanCell
# 
# 
# plot( habitat$polygon, col = "gray80", border = grey(0.3))
# plot( meanDensity.R, add = T,
#       breaks = cuts, col = col,
#       axes = F, box = F, bty = "n", legend = F)
# plot( habitat$polygon, add = T, border = grey(0.3))
# plot( st_geometry(countries), add = T, lwd = 2)
# 
# mtext( text = paste( "N = ", round(WA_SpaceUse$summary["Total",1],1),
#                      " [", round(WA_SpaceUse$summary["Total",4],1), " ; ",
#                      round(WA_SpaceUse$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# par(mfrow = c(1,1), mar = c(4,4,0,0))
# plot( Cleanpoly, col = "gray80", border = grey(0.3))
# plot( meanDensity.R, add = T,
#       breaks = cuts, col = col, colNA = NA,
#       axes = F, box = F, bty = "n", legend = F)
# plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
# mtext( text = paste( "N = ", round(SpaceUse$summary["Total",1],1),
#                      " [", round(SpaceUse$summary["Total",4],1), " ; ",
#                      round(SpaceUse$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# mtext( text = "Total utilization distribution", side = 3, font = 2, cex = 2, line = -2)
# 
# ##-- Add density legend
# plot( meanDensity.R, breaks = cuts, col = col,
#       legend.width = 1, legend.only = TRUE,
#       axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
#                         labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
#                         cex.axis = 0.6),
#       legend.args = list( text = 'Density (ids.km2)', line = -2,
#                           side = 4, font = 2, cex = 0.8))



##------------------------------------------------------------------------------