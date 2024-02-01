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



## ------   3. DNA DATA ------
ngs <- read.csv(file.path(dataDir,
                          "DNA/wa_genetic_samples_28012022_Pierre_SCR.csv"))
## Format dates and coordinates
ngs$Date <- as.POSIXct(strptime(ngs$Date, "%Y-%m-%d"))
ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
ngs$Month <- as.numeric(format(ngs$Date,"%m"))

ngs$x <- unlist(lapply(ngs$Coordinates.WGS84.UTM.zone.32N,
                       function(x){
                         as.numeric(strsplit(x, ", ")[[1]][1])
                       }))
ngs$y <- unlist(lapply(ngs$Coordinates.WGS84.UTM.zone.32N,
                       function(x){
                         as.numeric(strsplit(x, ", ")[[1]][2])
                       }))
## Turn into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$x,ngs$y)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"

## Plot check
plot(ngs, add = T, col = as.numeric(ngs$Temp.ID))

## number of detections per individual
table(ngs$Temp.ID, useNA = "always")

## number of detections per pack
table(ngs$Pack, useNA = "always")

## mean number of detections per sex
table(ngs$Sex, useNA = "always")
apply(table(ngs$Temp.ID, ngs$Sex, useNA = "always"), 2, function(x)sum(x)/sum(x>0))

## mean number of detections per status
table(ngs$Status, useNA = "always")
apply(table(ngs$Temp.ID,ngs$Status,useNA = "always"), 2, function(x)sum(x)/sum(x>0))




## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. HABITAT ------
## ------     1.1. HABITAT CHARACTERISTICS ------
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





## ------     1.2. HABITAT COVARIATES ------
## ------       1.2.1. CORINE LAND COVER ------
## Load data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))
plot(CLC)
plot(st_geometry(habitat$grid), add = T)

## Legend : 
## 1	= developed; 2 = agriculture; 3 = forest; 4 = herbaceous; 5 =	bare rock;
## 6 = sea features; 7 = inland water; 9 = perpetual snow
par(mfrow = c(1,2))
CLC.df <- matrix(NA, habitat$n.HabWindows, 9)
for(this.c in 1:9){
  temp <- CLC
  temp[!temp[] %in% this.c] <- 0
  temp[temp[] %in% this.c] <- 1
  plot(temp)
  plot(st_geometry(habitat$grid), add = T)
  
  COV <- raster::aggregate( x = temp,
                            fact = habitat$resolution/res(temp),
                            fun = mean)
  plot(COV)
  plot(st_geometry(habitat$grid), add = T)
  CLC.df[ ,this.c] <- raster::extract( x = COV,
                                       y = habitat$sp)
  print(this.c)
}#c


## Store corine land cover categories in a dataframe
CLC.df <- cbind.data.frame( "id" = 1:nrow(CLC.df),
                            "developed" = CLC.df[ ,1],
                            "agriculture" = CLC.df[ ,2],
                            "forest" = CLC.df[ ,3],
                            "herbaceous" = CLC.df[ ,4],
                            "bare rock" = CLC.df[ ,5],
                            "sea features" = CLC.df[ ,6],
                            "inland water" = CLC.df[ ,7],
                            "perpetual snow" = CLC.df[ ,9])

## Store in habitat grid
habitat$grid  <- habitat$grid %>% 
  left_join(CLC.df, by = "id") 




## ------       1.2.2. ELEVATION -------
## Load data
DEM <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))
plot(DEM)
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
DEM <- raster::aggregate( x = DEM,
                          fact = habitat$resolution/res(DEM),
                          fun = mean)
plot(DEM)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled elevation values
habitat$grid$elev <- DEM %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       1.2.3. TERRAIN RUGGEDNESS INDEX -------
## Load data
TRI <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))
plot(TRI)
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
TRI <- raster::aggregate( x = TRI,
                          fact = habitat$resolution/res(TRI),
                          fun = mean)
plot(TRI)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled terrain ruggedness 
habitat$grid$tri <- TRI %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       1.2.4. HUMAN POPULATION DENSITY -------
## Load data
POP <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
plot(POP)
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
POP <- raster::aggregate( x = POP,
                          fact = habitat$resolution/res(POP),
                          fun = mean)
plot(POP)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled human pop density 
habitat$grid$pop <- POP %>% 
  raster::extract(y = habitat$sp) %>% 
  scale()




## ------       1.2.5. ROAD DENSITY -------
## Load data
roads <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Road Density/road_density_tot_clipped.shp"))
table(roads$highway) # see https://wiki.openstreetmap.org/wiki/Key:highway for road classes description
mainRoads <- subset(roads, highway %in% c("motorway", "trunk", "primary"))
rm(roads)

## Extract road length in each habitat grid cell
intersection <- st_intersection(habitat$grid, mainRoads) %>%
  mutate(length = st_length(.))  %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(mainRoads_L = sum(length)) 

## Store scaled road density 
habitat$grid <- habitat$grid %>%
  left_join(intersection, by = "id") 
habitat$grid$mainRoads_L[is.na(habitat$grid$mainRoads_L)] <- 0
habitat$grid$mainRoads_L <- scale(habitat$grid$mainRoads_L)


## Display all habitat covariates
plot(habitat$grid, max.plot = 12)




## ------   2. DETECTORS ------
## For now: detector grid based on transect locations with the chosen resolution
## Snow fall, transect length an number of visits, is extracted for each detector cell
## each month. Only detector grid cells with counts > 0 are kept (i.e. where people logged GPS tracks)

## Alternative for faster model:
## Collapse all months and use average snowfall, total search length and total number of visits


## ------     2.1. DETECTORS CHARACTERISTICS ------
## Extract date and month for each transect
transects$Date <- as.numeric(str_sub(transects$Name, start= -6))
transects$Month <- as.numeric(str_sub(transects$Name, start= -4, end = -3))

## Create a grid from transects with the desired resolution
grid <- st_as_sf(st_make_grid(transects,detectors$detResolution))
grid$id <- 1:nrow(grid)

## Extract length and number of transects 
intersection <- st_intersection(grid, transects) %>%
  mutate(LEN = st_length(.))  %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
            transect_N = length(unique(Date)))   ## Get total number of visits in each detector grid cell
# head(intersection)
# plot(intersection$transect_L,intersection$transect_N)

## Store in detector grid
grid <- grid %>% 
  left_join(intersection, by = "id") %>% 
  filter(!is.na(transect_L))
plot(grid)
detectors$grid <- grid

## Extract detector coordinates
detectors$coords <- st_coordinates(st_centroid(grid))[ ,c("X","Y")]

dimnames(detectors$coords) <- list(1:nrow(detectors$detectors),
                                   c("x","y"))

##---- Extract number of detectors per month
detectors$n.detectors <- nrow(detectors$grid)



## ------     2.1. DETECTORS CHARACTERISTICS ------
# ## Create a list of detector grids (one for each month) 
# ## and extract grid cells that were searched each month with their associated 
# ## transect length and number of visits.
# for(m in 1:length(detectors$samplingMonths)){
#   ## Subset to transects searched in sampling month "m"
#   temp <- transects[transects$Month %in% detectors$samplingMonths[m], ]
#   
#   ## Create a grid from transects with the desired resolution
#   grid <- st_as_sf(st_make_grid(temp,detectors$detResolution))
#   grid$id <- 1:nrow(grid)
#   
#   ## Extract length and number of transects 
#   intersection <- st_intersection(grid, temp) %>%
#     mutate(LEN = st_length(.))  %>%
#     st_drop_geometry() %>%
#     group_by(id) %>%
#     summarise(transect_L = sum(LEN),
#               transect_N = length(LEN))
#   
#   ## Store in detector grid
#   grid <- grid %>% 
#     left_join(intersection, by = "id") %>% 
#     filter(!is.na(transect_L))
#   plot(grid)
#   detectors$grid.list[[m]] <- grid
#   
#   ## Extract detector coordinates
#   detectors$coords.list[[m]] <- st_coordinates(st_centroid(grid))[ ,c("X","Y")]
# }#m
# 
# 
# ##---- Extract number of detectors per month
# detectors$n.detectors <- unlist(lapply(detectors$grid.list,nrow))
# 
# ##---- Format detector coordinates in a rugged array
# detectors$coords <- array(NA, c(max(detectors$n.detectors),
#                                                   2,
#                                                   length(detectors$samplingMonths)))
# for(m in 1:length(detectors$samplingMonths)){
#   detectors$coords[1:detectors$n.detectors[m], ,m] <- detectors$coords.list[[m]]
# }#m
# 
# dimnames(detectors$coords) <- list(1:max(detectors$n.detectors),
#                                    c("x","y"),
#                                    detectors$samplingMonths)




## ------     2.2. DETECTOR COVARIATES ------
## ------       2.2.1. SNOWFALL ERA-5 ------
tif <- read_stars(file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/ERA-5_snowfall_alps_32N.nc"))
SNOW <- st_as_sf(tif)

for(m in 1:length(detectors$samplingMonths)){
  ## Subset to snow in month "m"
  temp <- SNOW[paste0("Band",m)]
  names(temp) <- c("cov","geometry")
  
  ## Extract snow fall in each detector grid cell
  intersection <- st_intersection(detectors$grid, temp) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(snow_fall = mean(cov)) 
  
  ## Store scaled road density 
  detectors$grid.list[[m]] <- detectors$grid.list[[m]] %>%
    left_join(intersection, by = "id") 
  detectors$grid.list[[m]]$snow_fall[is.na(detectors$grid.list[[m]]$snow_fall)] <- 0
  plot(detectors$grid.list[[m]])
}#m


# oct20 <- brick(f, var="Band1", lvar=4)
# nov20 <- brick(f, var="Band2", lvar=4)
# dec20 <- brick(f, var="Band3", lvar=4)
# jan21 <- brick(f, var="Band4", lvar=4)
# feb21 <- brick(f, var="Band5", lvar=4)
# mar21 <- brick(f, var="Band6", lvar=4)
# apr21 <- brick(f, var="Band7", lvar=4)
#
# plot(oct20)
# plot(transects,col = "blue", add=T)
# 
# ext_octsf <- raster::extract(snow, transects, df = T)
# mean_octsf <- ext_octsf %>%
#   setNames(c("ID", "sf_oct")) %>%          # rename for ease
#   group_by(ID) %>%                         # group by ID
#   summarise(mean = mean(sf_oct), n = n())  # calculate mean for cell (ID)



SNOW2 <-  st_drop_geometry(SNOW)
dim(SNOW2)

plot(unlist(SNOW2[1, ]), type = "l", col = 1, ylim = c(0,0.003))
for(i in seq(2, nrow(SNOW2),10)){
  points(unlist(SNOW2[i, ]), type = "l", col = i)
}




##---- Format detector coordinates in a rugged array
detectors$covs <- array(NA, c(max(detectors$n.detectors),
                              2,
                              length(detectors$samplingMonths)))
for(m in 1:length(detectors$samplingMonths)){
  detectors$covs[1:detectors$n.detectors[m], ,m] <- cbind(detectors$grid.list[[m]]$transect_L,
                                                          detectors$grid.list[[m]]$snow_fall) 
}#m



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
localObjects <- list()
for(m in 1:length(detectors$samplingMonths)){
  localObjects[[m]] <- getLocalObjects(
    habitatMask = habitat$binary,
    coords = scaledCoords$coordsDataScaled[ ,1:2,m],
    dmax = 21)
}#m


##---- Collapse number of local traps for each month into a matrix
localTrapsNum <- do.call(cbind,lapply(localObjects,function(x)x$numLocalIndices))
localTrapsNumMax <- unlist(lapply(localObjects,function(x)x$numLocalIndicesMax))


##---- Format local trap indices in a rugged array
localTrapsIndices <- array(NA, c(habitat$n.HabWindows,
                                 max(localTrapsNumMax),
                                 length(detectors$samplingMonths)))
for(m in 1:length(detectors$samplingMonths)){
  localTrapsIndices[ ,1:localTrapsNumMax[m],m] <- localObjects[[m]]$localIndices
}#m
dim(localTrapsIndices)




## ------   4. DETECTION DATA ------
## ------     4.1. DETECTION MATRIX : y ------
ngs$Id <- ngs$Temp.ID
ngs$Date
ngs$detectors <- apply(cbind(ngs$x,ngs$y), 1, function(x){
  getWindowIndex( curCoords = x,
                  lowerCoords = detectors$lowerCoords[,1:2],
                  upperCoords = detectors$upperCoords[,1:2]
  )
})

y.ar <- MakeY( myData = ngs,
               myDetectors = detectors$coords[,1:2],
               method = "Bernoulli",
               returnIdvector = TRUE)

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
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] ~ dnorm(0.0,0.01)
  }
  
  for(h in 1:n.habWindows){
    habIntensity[h] <- exp(inprod(betaHab[1:n.habCovs],
                                  hab.covs[h,1:n.habCovs]))
  }
  
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
  for(c in 1:n.detCovs){
    betaDet[c] ~ dnorm(0.0,0.01)
  }
  
  for(s in 1:n.states){
    for(ss in 1:2){
      p0[s,ss] ~ dunif(0,1)
      sigma[s,ss] ~ dunif(0,6)
      
      for(j in 1:n.detectors[k]){
        for(k in 1:n.occasions){
        logit(p0Traps[s,ss,j,k]) <- logit(p0[s,ss]) +
          inprod(betaDet[1:n.detCovs], det.covs[j,1:n.detCovs,k])
      }#j
    }#k
  }#ss
}#s


for(i in 1:n.individuals){
  for(k in 1:n.occasions){
    y[i,1:n.maxDets,k] ~ dbinomLocal_normal(
      size = size[1:n.detectors[k]],
      p0Traps = p0Traps[state[i],sex[i]+1,1:n.detectors[k],k],
      sigma = sigma[state[i],sex[i]+1],
      s = sxy[i,1:2],
      trapCoords = detCoords[1:n.detectors[k],1:2],
      localTrapsIndices = localTrapsIndices[1:n.habWindows,1:n.localIndicesMax],
      localTrapsNum = localTrapsNum[1:n.habWindows,k],
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      indicator = z[i],
      lengthYCombined = n.maxDets)
  }#k
}#i


##-- DERIVED PARAMETERS 
N <- sum(z[1:n.individuals])
})




## ------   2. DATA SIMULATION ------
##---- Set model constants, data & parameter simulated values (==inits)
simConstants <- list( n.individuals = 500,
                      n.habWindows = habitat$n.HabWindows,
                      n.habCovs = 4,
                      n.states = 3,
                      n.detectors = detectors$n.detectors, 
                      n.occasions = length(detectors$samplingMonths),
                      n.detCovs = 2,
                      n.localIndicesMax = localTrapsNumMax,
                      n.maxDets = 40,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

simData <- list( lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords, 
                 hab.covs = st_drop_geometry(habitat$grid[ ,c("forest","herbaceous","elev","mainRoads_L")]),
                 alpha = matrix(1,3,2),
                 det.covs = detectors$covs,
                 size = rep(1,max(detectors$n.detectors)),
                 detCoords = detectors$scaledCoords,
                 localTrapsIndices = localTrapsIndices,
                 localTrapsNum = localTrapsNum,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects[[1]]$habitatGrid)

simInits <- list( "betaHab" = c(0,0,-0.5,-0.2),
                  "psi" = 0.8,
                  "rho" = 0.5,
                  "theta" = matrix(c(0.5,0.45,0.05,
                                     0.5,0.3,0.2),3,2),
                  "betaDet" = c(0,0.05),
                  "p0" = matrix(c(0.2,0.1,0.05,
                                  0.2,0.1,0.05),3,2),
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
nodesToSim <- CsimModel$getDependencies(c(names(simInits)),
                                        self = F,
                                        includeData = F,
                                        downstream = T,
                                        returnScalarComponents = T)

nodesToSim.simple <- unique(sapply(strsplit(nodesToSim, "\\["), "[", 1))
nodesToSim.simple

##---- Simulate from NIMBLE 
CsimModel$simulate(nodesToSim)
CsimModel$calculate()



## ------   3. BUNDLE NIMBLE DATA/INITIAL VALUES/... ------
nimConstants <- simConstants 
nimData <- simData
nimInits <- simInits

##---- Identify which individuals were detected
whichDet <- rowSums(CsimModel$y[ ,1, ]) > 0
whichNotDet <- rowSums(CsimModel$y[ ,1, ]) == 0

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

state <- state.init <- CsimModel$sex
state[whichNotDet] <- NA
state.init[whichDet] <- NA
nimData$state <- state
nimInits$state <- state.init

nimInits$sxy <- CsimModel$sxy

nimParams <- c("N", "betaHab", "betaDet",
               "p0", "sigma", "psi", "rho", "theta")
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
  print(system.time(
    runMCMCbites( mcmc = Cmcmc,
                  bite.size = 500,
                  bite.number = 5,
                  path = file.path(simDir, modelName, paste0("output/chain",c))) 
  ))
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
