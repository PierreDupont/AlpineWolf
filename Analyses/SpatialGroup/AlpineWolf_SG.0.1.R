## -------------------------------------------------------------------------- ##
## ------------------------ ALPINE WOLF SG ---------------------------------- ##
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
modelName = "AlpineWolf.SG.0.1"
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
detectors <- ct
detectors$coords <- st_coordinates(detectors)
dimnames(detectors$coords) <- list(1:nrow(detectors$coords),
                                   c("x","y"))

detectors$`non-workin`[is.na(detectors$`non-workin`)] <- 0

n.detectors <- nrow(detectors)

## NEEDS TO FIGURE OUT HOW TO EXTRACT NUMBER OF DAYS, WEEKS AND MONTHS ACTIVE 
## FOR EACH CAMERA
detectors$months_active <- trunc(detectors$tot_attivi/30)+1
# detectors <- detectors %>% 
#   mutate(
#     "days_active" = lubridate::day(date_en) - lubridate::day(date_st),
#     "weeks_active" = lubridate::week(date_en) - lubridate::week(date_st),
#     "months_active" = lubridate::month(date_en-date_st)))
# 
# hist(detectors$days_active)
# hist(detectors$weeks_active)
hist(detectors$months_active)


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
##---- Calculate distance between dpictures and camera-traps
closest <- nn2( st_coordinates(detectors),
                st_coordinates(pics),
                k = 1,
                searchtype = "radius",
                radius = 30000)
table(closest$nn.dists)

##---- Identify and plot pictures to check
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
##---- Aggregate pictures per DAY
pics_d <- pics %>%
  group_by(day = lubridate::day(date), CT_id) %>%
  summarise(n_wolves = max(n_wolves))
hist(pics_d$n_wolves)
dim(pics_d)

##---- Aggregate pictures per WEEK
pics_w <- pics %>%
  group_by(week = lubridate::week(date), CT_id) %>%
  summarise(n_wolves = max(n_wolves))
hist(pics_w$n_wolves)
dim(pics_w)

##---- Aggregate pictures per MONTH
pics_m <- pics %>%
  group_by(Month, CT_id) %>%
  summarise(n_wolves = max(n_wolves))
hist(pics_m$n_wolves)
dim(pics_m)

par(mfrow=c(2,2))
hist(pics$n_wolves, breaks = seq(0.5,12.5,by = 1))
hist(pics_w$n_wolves, breaks = seq(0.5,12.5,by = 1))
hist(pics_w$n_wolves, breaks = seq(0.5,12.5,by = 1))
hist(pics_m$n_wolves, breaks = seq(0.5,12.5,by = 1))




## ------     5.3. DETECTION VECTOR : y ------
months <- sort(na.omit(unique(pics$Month)))
n.months <- length(months)

##-- Create a vector of camera-trap specific observations
y <- y.visits <- matrix(0, nrow = n.detectors, ncol = n.months)

##-- Loop over camera-traps
for(j in 1:nrow(detectors)){
  ##-- Identify this detector
  thisCT <- detectors$id[j]
  
  theseMonths <- month(detectors$date_st[i])
  
  
  for(m in 1:n.months){
    
    ##-- Subset pictures to this camera-trap
    thesePics <- pics %>%
      filter(.,
             CT_id == thisCT, 
             Month == months[m])
    
    if(nrow(thesePics)>0){
      ##-- Calculate total number of wolves detected at this camera-trap
      y[j] <- max(thesePics$n_wolves)
      ##-- Calculate total number of pictures (=events) at this camera-trap
      y.visits[j] <- nrow(thesePics)
    }#if
  }#i
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
  
  ##-- Loop over groups
  for(g in 1:G) {
    s[g,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##---- DEMOGRAPHIC PROCESS  
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
    for(k in 1:K[j]){
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



## ------   2. BUNDLE DATA ------
##---- Set model data
nimData <- list( y = y,
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
nimConstants <- list( G = 500,
                      J = nrow(y),
                      
                      n.habWindows = habitat$n.HabWindows,
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2],
                      n.habCovs = ncol(nimData$hab.covs))

##---- Set parameters to track 
nimParams <- c("N", "D", "lambda0", "sigma", "psi", "theta", "rho", "betaHab")
nimParams2 <- c("z", "s")




## ------   3. NIMBLE INITIAL VALUES ------
##---- Generate initial individual AC locations
sInits <- matrix(NA, nrow = M, ncol = 2)
for (i in 1:M) {
  sInits[i, ] <- rbernppAC( 
    n = 1,
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
    sigma = matrix(c(0.7334833, 0.4788753,
                     0.5919603, 0.5573514,
                     2.6903494, 5.9257740), byrow = T, nrow = 3),
    alpha0 = 3,
    lambda0 = matrix(c( 0.88,2.08,
                        1.36, 1.54,
                        0.07,0.02), byrow = T, nrow = 3),
    oper = operInits,
    s = sInits,
    lambda.oper = 1,
    z = rbinom(M,1,0.6),
    rho = 0.5,
    theta = matrix(c(0.2784087, 0.29532485,
                     0.6189801, 0.61738292,
                     0.1026112, 0.08729223), byrow = T, nrow = 3),
    sex = rbinom(M,1,0.5),
    status = rcat(n = M, prob = c(0.28,0.62,0.1)),
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



