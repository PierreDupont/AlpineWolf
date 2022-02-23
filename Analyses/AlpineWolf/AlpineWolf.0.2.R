################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[.]
##### ---------------- z[psi,rho,theta]
##### ---------------- y[p0,sigma] ---------------------- #####
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
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity_PD.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.0.2"
thisDir <- file.path(analysisDir, modelName)

## HABITAT SPECIFICATIONS
habitat = list( country =  c("SWE","NOR"),
                resolution = 5000, 
                buffer = 20000)

## DETECTORS SPECIFICATIONS
detectors = list( detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

## NGS DATA SPECIFICATIONS
data = list( sex = c("F","M"),
             status = c("alpha","pup","other"),
             aug.factor = 6) 

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}



## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
## NEED COUNTRY LEVEL POLYGON TO CUT OUT SEA AND OTHER NON-HABITAT
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() 

## Plot check
plot(st_geometry(studyArea))
plot(st_geometry(countries), col = "gray80",add=T)
plot(studyAreaGrid["SCR"], add = T)



## ------   2. SEARCH EFFORT DATA ------
## Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))

## Convert dates
transects$Date <- parse_date_time(transects$date, orders = c('ymd'))
transects$Year <- as.numeric(format(transects$Date,"%Y"))
transects$Month <- as.numeric(format(transects$Date,"%m"))

## Plot check
plot(transects, col = "red", add = T)




## ------   3. DNA DATA ------
# allSamples <- read_sf(file.path(dataDir,"GISData/scats/merged_scats.shp"))
# plot(st_geometry(studyArea))
# plot(st_geometry(countries), col = "gray80",add=T)
# plot(allSamples, add = T,pch=3)

ngs <- read.csv(file.path(dataDir,"DNA/ngs.csv"))
dim(ngs)

## Use lubridate to clean dates
# The orders argument is a character vector containing the possible date-time parsing 
# formats in the order they should be tested. So by giving c('dmy', 'ymd'), lubridate 
# will try to parse all strings as Day, Month, Year format. If it can't do that 
# successfully (for example, the date 2021-03-18 won't work as there is no 2021th day), 
# it will try the next in the list until all strings are parsed, or all orders exhausted.
ngs$Date <- parse_date_time(ngs$Date, orders = c('dmy','mdy','ymd'))
ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
ngs$Month <- as.numeric(format(ngs$Date,"%m"))

## Filter out samples from 2019
ngs <- ngs[ngs$Year > 2019, ]
dim(ngs)

## Filter out dead recoveries
ngs <- ngs[ngs$Dead.recovery != "Yes", ]
dim(ngs)

## Number of detections per individual
numDetsPerId <- table(ngs$Genotype.ID)
hist(numDetsPerId)

## Number of individuals detected
length(numDetsPerId)

## Mean number of detections per ID
mean(numDetsPerId)

## Number of individuals with > 1 detection
sum(numDetsPerId > 1)

## mean number of individual detections by sex
table(ngs$Sex, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Sex, useNA = "always"),
      2,
      function(x)sum(x)/sum(x>0))

## mean number of individual detections by status
ngs$Status[ngs$Status == "pup 2016"] <- "pup"
ngs$Status[ngs$Status == "pup 2017"] <- "pup"
table(ngs$Status, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Status, useNA = "always"),
      2,
      function(x)sum(x)/sum(x>0))

## mean number of individual detections by sex and status
apply(table(ngs$Genotype.ID, ngs$Status, ngs$Sex, useNA = "always"),
      c(2,3),
      function(x)sum(x)/sum(x>0))

## Turn ngs into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$Coordinate.WGS84.UTM.East,
                                     ngs$Coordinate.WGS84.UTM.North)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

## Plot check
plot(st_geometry(studyArea))
plot(st_geometry(countries), col = "gray80",add=T)
plot(st_geometry(transects), col = "red", add = T)
plot(st_geometry(ngs), add = T, pch = 3)

# ## Calculate distance to the closest transect
# ## For later: can be quite long with many points and transects
# dist <- st_distance(ngs,transects)
# ngs$dist_to_transect <- apply(dist,1,min)
# range(ngs$dist_to_transect)

## Fix one problematic wolf
ngs$Genotype.ID[ngs$Genotype.ID=="WBS-M002"][5:7] <- "WBS-M002.2"
ngs[ngs$Genotype.ID=="WBS-M002", ]


## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. DETECTORS ------
## ------     1.1. DETECTORS CHARACTERISTICS ------
# ## Create a line to delineate the "western Alps"
# cutline <- st_linestring(rbind(c(630000,4800000),
#                                c(430000,5230000)))
# cutline <- st_sfc(cutline)
# st_crs(cutline) <- st_crs(studyArea)
# cutline <- st_buffer(cutline, dist = 0.0001)
# plot(studyArea)
# plot(cutline, add = T, border = "red", lwd = 2)
# 
# ##---- Create study area polygon for the western Alps (study area cut to the line) 
# studyArea_west <- st_sfc(st_cast(st_difference(st_buffer(studyArea,500), cutline),"POLYGON")[[1]])
# st_crs(studyArea_west) <- st_crs(studyArea)
# studyArea_west <- st_intersection( st_union(countries),
#                                    studyArea_west,
#                                    drop_lower_td = TRUE)
# plot(studyArea_west, add=T, col = "gray60")

##---- Create detector raster 
detectors$raster <- raster(extent(st_bbox(studyArea)))
res(detectors$raster) <- detectors$detResolution
detectors$raster <- fasterize( countries, detectors$raster)

##---- Mask and crop raster cells outside the study area to obtain the detector grid
detectors$raster <- mask(detectors$raster, st_as_sf(studyArea))
detectors$raster <- crop(detectors$raster, st_as_sf(studyArea))
detectors$grid <- st_as_sf(rasterToPolygons(detectors$raster))
detectors$grid$id <- 1:nrow(detectors$grid)
st_crs(detectors$grid) <- st_crs(studyArea)

##---- Extract length and number of transects in each grid cell
intersection <- st_intersection(detectors$grid, transects) %>%
  mutate(LEN = st_length(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
            transect_N = length(unique(Date)))   ## Get total number of visits in each detector grid cell

##---- Store in detector grid
detectors$grid <- detectors$grid %>% 
  left_join(intersection, by = "id") 
detectorsgrid$transect_L[is.na(grid$transect_L)] <- 0
detectorsgrid$transect_N[is.na(grid$transect_N)] <- 1


####---- If you want only grid cells that contain transects instead
detectors$grid <- detectors$grid %>% 
  left_join(intersection, by = "id") %>% 
 filter(!is.na(transect_L))

detectors$grid$mean_transect_L <- scale(detectors$grid$transect_L/detectors$grid$transect_N)
detectors$trials <- detectors$grid$transect_N

##---- Extract detector coordinates
detectors$coords <- st_coordinates(st_centroid(detectors$grid))[ ,c("X","Y")]
dimnames(detectors$coords) <- list(1:nrow(detectors$grid),
                                   c("x","y"))
detectors$sp <- SpatialPoints(coords = detectors$coords,
                            proj4string = crs(detectors$grid))


##---- Extract total number of detectors
detectors$n.detectors <- nrow(detectors$grid)




## ------     1.2. DETECTOR COVARIATES ------
## ------       1.2.1. SNOWFALL ERA-5 ------
tif <- read_stars(file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/ERA-5_snowfall_alps_32N.nc"))
SNOW <- st_as_sf(tif)
SNOW$snow.mean <- rowMeans(st_drop_geometry(SNOW), na.rm = T) 
SNOW$snow.sum <- rowSums(st_drop_geometry(SNOW), na.rm = T) 

## Extract snow fall in each detector grid cell
intersection <- st_intersection(detectors$grid, SNOW) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(snow_fall = mean(snow.mean)) 

## Store scaled road density 
detectors$grid <- detectors$grid %>%
  left_join(intersection, by = "id") 

detectors$grid$snow_fall[is.na(detectors$grid$snow_fall)] <- 0
detectors$grid$snow_fall <- scale(detectors$grid$snow_fall)




## ------       1.2.2. EAST/WEST ------
## Create a line to delineate the "western Alps"
cutline <- st_linestring(rbind(c(630000,4800000),
                               c(430000,5230000)))
cutline <- st_sfc(cutline)
st_crs(cutline) <- st_crs(studyArea)
cutline <- st_buffer(cutline, dist = 0.0001)
plot(studyArea)
plot(cutline, add = T, border = "red", lwd = 2)

##---- Create study area polygon for the western Alps (study area cut to the line) 
studyArea_split <- st_difference(st_buffer(studyArea,1000), cutline)
studyArea_west <- st_intersection( st_union(countries),
                                   studyArea_west,
                                   drop_lower_td = TRUE)

intersection <- st_intersection(detectors$grid, studyArea_split) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(snow_fall = mean(snow.mean)) 

plot(studyArea_west, add=T, col = "gray60")

## ------       1.2.1. CORINE LAND COVER ------
## Load data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))
plot(CLC)
plot(st_geometry(detectors$grid), add = T)

## Legend : 
## 1	= developed; 2 = agriculture; 3 = forest; 4 = herbaceous; 5 =	bare rock;
## 6 = sea features; 7 = inland water; 9 = perpetual snow
layer.names <- c("developed","agriculture","forest","herbaceous","bare rock",
                 "sea features","inland water","NA","perpetual snow")
layer.index <- c(1,2,3,4,5,9)
par(mfrow = c(1,2))

CLC.df <- as.data.frame(matrix(NA, detectors$n.detectors, length(layer.index)))
names(CLC.df) <- layer.names[layer.index]
CLC.df$id <- 1:detectors$n.detectors

COV <- list()
for(l in 1:length(layer.index)){
  temp <- CLC
  temp[!temp[] %in% layer.index[l]] <- 0
  temp[temp[] %in% layer.index[l]] <- 1
  plot(temp)
  plot(st_geometry(detector$grid), add = T)
  
  COV[[l]] <- raster::aggregate( x = temp,
                                 fact = detector$resolution/res(temp),
                                 fun = mean)
  plot(COV[[l]])
  plot(st_geometry(detector$grid), add = T)
  CLC.df[ ,l] <- raster::extract( x = COV[[l]],
                                  y = detectors$sp)
  print(layer.index[l])
}#c

## Store in habitat grid
detectors$grid <- left_join(detectors$grid, CLC.df, by = "id") 




## ------       1.2.2. ELEVATION -------
## Load data
DEM_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))
plot(DEM_raw)
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
DEM <- raster::aggregate( x = DEM_raw,
                          fact = habitat$resolution/res(DEM_raw),
                          fun = mean)
plot(DEM)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled elevation values
habitat$grid$elev <- DEM %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       1.2.3. TERRAIN RUGGEDNESS INDEX -------
## Load data
TRI_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))
plot(TRI_raw )
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
TRI <- raster::aggregate( x = TRI_raw ,
                          fact = habitat$resolution/res(TRI_raw ),
                          fun = mean)
TRI <- raster::focal(x = TRI,
                     w = matrix(1,3,3),
                     fun = mean,
                     na.rm = T)
plot(TRI)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled terrain ruggedness 
habitat$grid$tri <- TRI %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       1.2.4. HUMAN POPULATION DENSITY -------
## Load data
POP_raw  <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
plot(POP_raw )
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
POP <- raster::aggregate( x = POP_raw ,
                          fact = habitat$resolution/res(POP_raw ),
                          fun = mean)
POP <- raster::focal(x = POP,
                     w = matrix(1,3,3),
                     fun = mean,
                     na.rm = T)
plot(POP)
plot(st_geometry(habitat$grid), add = T)

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



## ------       2.2.5. ROAD DENSITY 


## ------   2. HABITAT ------
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (study area + buffer)  
## Buffer the detector grid to delineate the habitat polygon
habitat$polygon <- st_buffer(st_union(detectors$grid),
                             habitat$buffer)

## Remove un-suitable habitat (e.g. seas)
habitat$polygon <- st_intersection( st_union(countries),
                                    habitat$polygon,
                                    drop_lower_td = TRUE)

##---- Create "ROOT" habitat raster (keep raster cells with > 50% habitat) 
habitat$raster <- raster(extent(st_bbox(st_buffer(st_union(studyAreaGrid),
                                                  100000))))
## Trick to get % of habitat
res(habitat$raster) <- habitat$resolution/10
habitat$raster <- fasterize( countries,
                             habitat$raster)
habitat$raster[is.na(habitat$raster)]<- 0
habitat$raster <- aggregate(habitat$raster, fact = 10, sum, na.rm = TRUE)
habitat$raster[habitat$raster < 50] <- NA
habitat$raster[habitat$raster >= 50] <- 1

# mask and crop habitat with the buffer 
habitat$raster <- mask(habitat$raster, st_as_sf(habitat$polygon))
habitat$raster <- crop(habitat$raster, st_as_sf(habitat$polygon))

habitat$grid <- st_as_sf(rasterToPolygons(habitat$raster,
                                          fun = function(x)x == 1))
habitat$grid$id <- 1:nrow(habitat$grid)
habitat$grid <- habitat$grid[,c(2,3)]
st_crs(habitat$grid) <- st_crs(habitat$polygon)

##---- Create "COMPARISON" habitat raster (for comparison between models) 
comparison <- st_buffer(st_union(detectors$grid[as.numeric(detectors$grid$transect_L) >0, ]),5000) 
habitat$comparison <- mask(habitat$raster, st_as_sf(comparison))
habitat$comparison[habitat$comparison == 0] <- NA

##---- Create "ITALY" habitat raster (for density estimates) 
Italia <- countries[countries$CNTR_CODE == "IT", ]
habitat$Italia <- mask(habitat$raster, st_as_sf(Italia))
habitat$Italia[habitat$Italia == 0] <- NA

##---- Create Habitat matrix of cell ID 
habitat$matrix <- habitat$raster
habitat$matrix[] <- 0
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
habitat$lowerCoords <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)

##---- Visual plotting to check if everything is right  
plot(habitat$raster)					
plot(st_geometry(countries), add = T)
#plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
#plot(studyArea_west, add = T, col = "gray40")
plot(transects,add=T, col="red")




## ------     2.2. HABITAT COVARIATES ------
## ------       2.2.1. CORINE LAND COVER ------
## Load data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))
plot(CLC)
plot(st_geometry(habitat$grid), add = T)

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




## ------       2.2.2. ELEVATION -------
## Load data
DEM_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))
plot(DEM_raw)
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
DEM <- raster::aggregate( x = DEM_raw,
                          fact = habitat$resolution/res(DEM_raw),
                          fun = mean)
plot(DEM)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled elevation values
habitat$grid$elev <- DEM %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       2.2.3. TERRAIN RUGGEDNESS INDEX -------
## Load data
TRI_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))
plot(TRI_raw )
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
TRI <- raster::aggregate( x = TRI_raw ,
                          fact = habitat$resolution/res(TRI_raw ),
                          fun = mean)
TRI <- raster::focal(x = TRI,
                     w = matrix(1,3,3),
                     fun = mean,
                     na.rm = T)
plot(TRI)
plot(st_geometry(habitat$grid), add = T)

## Extract scaled terrain ruggedness 
habitat$grid$tri <- TRI %>%
  raster::extract(y = habitat$sp) %>%
  scale()




## ------       2.2.4. HUMAN POPULATION DENSITY -------
## Load data
POP_raw  <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
plot(POP_raw )
plot(st_geometry(habitat$grid), add = T)

## Aggregate to the detector resolution
POP <- raster::aggregate( x = POP_raw ,
                          fact = habitat$resolution/res(POP_raw ),
                          fun = mean)
POP <- raster::focal(x = POP,
                     w = matrix(1,3,3),
                     fun = mean,
                     na.rm = T)
plot(POP)
plot(st_geometry(habitat$grid), add = T)

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



## ------       2.2.5. ROAD DENSITY -------
## Load data
roads  <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Road Density/road_density_tot_clipped.shp"))
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
  dmax = 9)



## ------   4. DETECTION DATA ------
## ------     4.1. DETECTION MATRIX : y ------
## Assign detector based on minimum distance
closest <- nn2( st_coordinates(st_centroid(detectors$grid)),
                st_coordinates(ngs),
                k = 1,
                searchtype = "radius",
                radius = 10000)
ngs$Detector <- c(closest$nn.idx)

# Count individual detections per detector
detMat <- as.matrix(table(ngs$Genotype.ID,ngs$Detector))

# Convert to binary detections
detMat[detMat > 0] <- 1

# RETRIEVE THE NUMBER OF DETECTIONS FOR EACH ID
detNums <- apply(detMat, 1, function(x) sum(x>0))

# INITIALIAZE EMPTY ARRAYS FOR THE DETECTIONS AND DETECTOR INDICES
detIndices <- matrix(-1, dim(detMat)[1], max(detNums)*2)
ySparse <- matrix(-1, dim(detMat)[1], max(detNums)*2)

# FILL IN THE ARRAYS
n.detected <- dim(detMat)[1]
for(i in 1:n.detected){
  # GET WHERE (DETECTOR ID) DETECTIONS OCCUR
  detIndices[i,1:detNums[i]] <- as.numeric(names(which(detMat[i, ] > 0)))
  # GET NUMBER OF DETECTIONS 
  ySparse[i,1:detNums[i]] <- detMat[i,which(detMat[i, ] > 0)]
}
yCombined <- cbind(detNums, ySparse, detIndices)




## ------     4.2. INDIVIDUAL COVARIATES ------
IDs <- dimnames(yCombined)[[1]]
sex <- status <- pack <- rep(NA, length(IDs))

for(i in 1:length(IDs)){
  ## Sex 
  temp <- unique(ngs$Sex[ngs$Genotype.ID %in% IDs[i]])
  sex[i] <- ifelse(length(temp>1),temp[1],temp)
  
  ## Social status
  temp <- unique(ngs$Status[ngs$Genotype.ID %in% IDs[i]])
  status[i] <-  ifelse(length(temp>1),temp[1],temp)
  
  # ## Pack membership
  # pack[i] <- unique(ngs$Pack[ngs$Genotype.ID %in% IDs[i]])
}

sex <- as.numeric(as.factor(sex)) - 1

status[status == "alpha"] <- 1
status[status == "pup"] <- 2
status[status == "other"] <- 3
status[status == ""] <- NA
status <- as.numeric(status)




## ------     4.3. DATA AUGMENTATION ------
yCombined.aug <- MakeAugmentation( y = yCombined,
                                   aug.factor = data$aug.factor,
                                   replace.value = 0)

sex.aug <- MakeAugmentation( y = sex,
                             aug.factor = data$aug.factor,
                             replace.value = NA)

status.aug <- MakeAugmentation( y = status,
                                aug.factor = data$aug.factor,
                                replace.value = NA)




## -----------------------------------------------------------------------------
## ------ III. VISUAL SUMMARY OF DATA PROCESSING ------
pdf( file = file.path(thisDir, paste0(modelName, "_data_prep.pdf")),
     width = 20, height = 12)


## ------   1. BEFORE/AFTER PLOTS ------
par(mfrow = c(1,2))

##-- RAW DATA
##-- Study area polygon 
plot(studyArea, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
mtext( text = "Raw data", side = 1, font = 2)
##-- GPS search tracks 
plot(transects, col = "firebrick3", add = T, lwd = 0.5)
##-- NGS samples per sex
sexCol <- ifelse(ngs$Sex == "F", "lightgreen", "orange4")
plot(ngs, add = T, pch = 21, bg = sexCol, col = sexCol, cex = 0.5)
legend( x = 372000, y = 5210000, bty = "n",
        legend = c("Females", "males"), title = "NGS samples",
        pch = 21, pt.bg = c("lightgreen", "orange4"))


##-- PROCESSED DATA
plot(habitat$polygon, border = "red", lwd = 2)
mtext( text = "Processed data", side = 1, font = 2)
##-- Habitat raster
temp <- habitat$raster
temp[temp[] == 0] <- NA
plot(temp,add=T, col = adjustcolor("forestgreen", alpha.f = 0.3), legend = F)
plot(st_geometry(countries),add=T)
plot(habitat$polygon, add = T, border = "red", lwd = 2)

##-- Main detectors
plot(st_centroid(detectors$grid), add = T, cex = 0.3, col = "gray40", pch = 3)
##-- Number of detections per detector
numDetsPerDet <- apply(detMat, 2, function(x)sum(x > 0))
colPal <- hcl.colors(n = max(numDetsPerDet)+1)
plot(st_centroid(detectors$grid[as.numeric(names(numDetsPerDet)), ]),
     add = T, pch = 21,
     cex = numDetsPerDet/2,
     bg = colPal[numDetsPerDet+1],
     col = colPal[numDetsPerDet+1])
mtext( text = "Data Processing", side = 3, line = -2, outer = TRUE, font = 2, cex = 2)

nameNumDets <- as.numeric(names(c(table(numDetsPerDet))))
legend( x = 502000, y = 5110000, bty = "n",
        legend = rev(nameNumDets), title = "# of obs. per detector",
        y.intersp = 1.2, x.intersp = 2,
        pch = 21, pt.bg = rev(colPal[nameNumDets+1]), pt.cex = rev((nameNumDets+1)/2))




## ------   2. DENSITY COVARIATES ------
par(mfrow = c(1,2), mar = c(2,2,2,2))
## ------     2.1. CORINE LAND COVER 
for(l in 1:length(layer.index)){
  temp <- CLC
  temp[!temp[] %in% layer.index[l]] <- 0
  temp[temp[] %in% layer.index[l]] <- 1
  
  plot(st_geometry(habitat$grid))
  mtext( text = paste("Raw",layer.names[layer.index[l]]),
         side = 1, font = 2)
  plot(temp, add = T)
  plot(st_geometry(countries), add = T)
  plot(st_geometry(habitat$grid), add = T)
  
  plot(st_geometry(habitat$grid))
  plot(habitat$grid[ ,layer.names[layer.index[l]]], add=T)
  mtext( text = paste("Processed",layer.names[layer.index[l]]),
         side = 1, font = 2)
  plot(st_geometry(countries), add = T)
  
  
  mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
}#c



## ------     2.2. ELEVATION 
## Raw data
plot(st_geometry(habitat$grid))
mtext( text = "Raw elevation",
       side = 1, font = 2)
plot(DEM_raw, add = T)
plot(st_geometry(countries), add = T)
plot(st_geometry(habitat$grid), add = T)

## Processed data
plot(st_geometry(habitat$grid))
mtext( text = "Processed elevation",
       side = 1, font = 2)
plot(habitat$grid[ ,"elev"], add = T)
plot(st_geometry(countries), add = T)

mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)



## ------     2.3. TERRAIN RUGGEDNESS INDEX 
## Raw data
plot(st_geometry(habitat$grid))
mtext( text = "Raw terrain ruggedness",
       side = 1, font = 2)
plot(TRI_raw, add = T)
plot(st_geometry(countries), add = T)
plot(st_geometry(habitat$grid), add = T)

## Processed data
plot(st_geometry(habitat$grid))
mtext( text = "Processed terrain ruggedness",
       side = 1, font = 2)
plot(habitat$grid[ ,"tri"], add = T)
plot(st_geometry(countries), add = T)

mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)


## ------     2.4. HUMAN POPULATION DENSITY 
par(mfrow = c(1,3), mar = c(2,2,2,2))

## Raw data
plot(st_geometry(habitat$grid))
mtext( text = "Raw human density",
       side = 1, font = 2)
plot(POP_raw, add = T)
plot(st_geometry(countries), add = T)
plot(st_geometry(habitat$grid), add = T)

## Processed data
plot(st_geometry(habitat$grid))
mtext( text = "Processed human density",
       side = 1, font = 2)
plot(habitat$grid[ ,"pop"], add = T)
plot(st_geometry(countries), add = T)

mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)

## Processed data
plot(st_geometry(habitat$grid))
mtext( text = "log(human density+1)",
       side = 1, font = 2)
plot(habitat$grid[ ,"log_pop"], add = T)
plot(st_geometry(countries), add = T)

mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)



## ------     2.5. ROAD DENSITY 
par(mfrow = c(1,2), mar = c(2,2,2,2))

## Raw data
plot(st_geometry(habitat$grid))
mtext( text = "Raw main roads",
       side = 1, font = 2)
plot(mainRoads[,"highway"], add = T)
plot(st_geometry(countries), add = T)
plot(st_geometry(habitat$grid), add = T)

## Processed data
plot(st_geometry(habitat$grid[ ,"mainRoads_L"]))
mtext( text = "Processed main roads",
       side = 1, font = 2)
plot(habitat$grid[ ,"mainRoads_L"], add = T)
plot(st_geometry(countries), add = T)

mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)



## ------   3. DETECTOR COVARIATES ------
par(mfrow = c(1,2), mar = c(2,2,2,2))

##-- Search transects
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
plot( transects, col = "firebrick3", add = T)
plot(st_geometry(countries), add = T)
mtext( text = "Search transects (raw)", side = 1, font = 2)

plot( habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
plot( detectors$grid[ ,"transect_L"], add = T)
plot(st_geometry(countries), add = T)
mtext( text = "Search transects length (processed)", side = 1, font = 2)

mtext( text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)


##-- Snow
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
plot( SNOW[ ,"snow.sum"], add = T)
plot(habitat$polygon, add=T)
plot(st_geometry(countries), add = T)
mtext( text = "Cumulative snow fall (raw)", side = 1, font = 2)

plot( habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
plot( detectors$grid[,"snow_fall"], add = T)
plot(st_geometry(countries), add = T)
mtext( text = "Cumulative snow fall(processed)", side = 1, font = 2)

mtext( text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)



## ------   4. INDIVIDUAL COVARIATES -----
##-- individual centroids
par(mfrow = c(1,1))
centroids <- st_drop_geometry(ngs) %>%
  group_by(Genotype.ID) %>%
  summarise(X = mean(Coordinate.WGS84.UTM.East),               ## Get total length searched in each detector grid cell
            Y =  mean(Coordinate.WGS84.UTM.North)) 
coordinates(centroids) <- cbind.data.frame(centroids$X,
                                           centroids$Y)
proj4string(centroids) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
centroids <- st_as_sf(centroids)
plot(habitat$polygon, col = "gray20", border = F)
plot(transects, col = "gray80", add = T)
plot(centroids,
     add = T, pch = 21, cex = 1,
     bg = adjustcolor("red",0.2),
     col = "red")

mtext( text = "Individual centroids", side = 3, font = 2)


##-- NGS samples per sex
par(mfrow = c(1,2))
sexCol <- c(hcl.colors(2),"gray60")
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
plot(transects, col = "black", add = T)
plot(st_geometry(ngs),
     add = T, pch = 21, cex = 1,
     bg = ifelse(ngs$Sex == "F", sexCol[1], sexCol[2]),
     col = "black")
legend( x = 502000, y = 5110000, bty = "n",
        legend = c("Females", "males", "NA"),
        title = "NGS samples",
        pch = 21,
        pt.bg = sexCol)

barplot( cbind(table(sex,useNA = "always")),
         col = sexCol,
         ylab = "Number of individuals detected",
         beside = F, width = 1, xlim = c(0,2))
legend( x = 1.3, y = 320,
        legend = c("NA", "male","female"),
        fill =  rev(sexCol))

##-- NGS samples per sex
statusCol <- c(hcl.colors(3),"gray60")
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
plot(transects, col = "black", add = T, lwd = 0.5)
plot(st_geometry(ngs),
     add = T, pch = 21, cex = 0.5,
     bg = ifelse(ngs$Status == "alpha", statusCol[1],
                 ifelse(ngs$Status == "pup",statusCol[2],
                        ifelse(ngs$Status == "other",statusCol[3],
                               statusCol[4]))),
     col = "black")
legend( x = 502000, y = 5110000, bty = "n",
        legend = c("alpha", "pup", "other", "NA"),
        title = "NGS samples",
        pch = 21, pt.bg = statusCol)

barplot( cbind(table(status,useNA = "always")),
         col = statusCol,
         beside = F, width = 1, xlim = c(0,2))
legend( x = 1.3, y =320,
        legend = c("NA", "other", "pup", "alpha"),
        fill =  rev(statusCol))

graphics.off()



## -----------------------------------------------------------------------------
## ------ IV. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] ~ dnorm(0.0,0.01)
  }#c
  for(h in 1:n.habWindows){
    habIntensity[h] <- exp(inprod(betaHab[1:n.habCovs],
                                  hab.covs[h,1:n.habCovs]))
  }#h
  
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
    betaDet ~ dnorm(0.0,0.01)
  }
  
  for(s in 1:n.states){
    for(ss in 1:2){
      p0[s,ss] ~ dunif(0,1)
      sigma[s,ss] ~ dunif(0,6)
      for(j in 1:n.detectors){
        logit(p0Traps[s,ss,j]) <- logit(p0[s,ss]) + inprod(betaDet[1:n.detCovs],
                                                           det.covs[h,1:n.detCovs])
      }#j
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




## ------   2. BUNDLE DATA ------
##---- Set model constants, data & parameter simulated values (==inits)
nimConstants <- list( n.individuals = dim(yCombined.aug)[1],
                      n.maxDets = dim(yCombined.aug)[2],
                      n.habWindows = habitat$n.HabWindows,
                      # habIntensity = rep(1,habitat$n.HabWindows),
                      n.detectors = detectors$n.detectors, 
                      n.habCovs = 6,
                      n.detCovs = 2,
                      n.states = 3,
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

nimData <- list( y = yCombined.aug,
                 z = c(rep(1,n.detected),
                       rep(NA,dim(yCombined.aug)[1]-n.detected)),
                 sex = sex.aug,
                 status = status.aug,
                 alpha = matrix(1,3,2),
                 lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords, 
                 hab.covs = cbind(habitat$grid$`log_pop`,
                                  habitat$grid$`roads`,
                                  habitat$grid$`forest`,
                                  habitat$grid$`herbaceous`,
                                  habitat$grid$`bare rock`,
                                  habitat$grid$`perpetual snow`),
                 det.covs = cbind(detectors$grid$`transect_L`,
                                  detectors$grid$`snow_fall`),
                 size = rep(1,detectors$n.detectors),# detectors$grid$transect_N,
                 detCoords = detectors$scaledCoords,
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects$habitatGrid)

nimParams <- c("N", "p0", "sigma", "psi",
               "betaDet", "betaHab", "theta", "rho",
               "z", "s")


## ------   3. SAVE THE INPUT ------
for(c in 1:4){
  s.init <- matrix(NA, nimConstants$n.individuals, 2)
  for(i in 1:n.detected){
    if(detNums[i] > 1){
      s.init[i, ] <- colMeans(detectors$scaledCoords[detIndices[i,1:detNums[i]], ])
    } else {
      s.init[i, ] <- detectors$scaledCoords[detIndices[i,1:detNums[i]], ] + rnorm(2,0,0.1)
    }
  }#i
  for(i in (n.detected + 1):nimConstants$n.individuals){
    s.init[i, ] <- rbernppAC( n = 1,
                              lowerCoords = nimData$lowerHabCoords,
                              upperCoords = nimData$upperHabCoords,
                              logIntensities = log(rep(1,habitat$n.HabWindows)),
                              logSumIntensity = log(sum(rep(1,habitat$n.HabWindows))),
                              habitatGrid = nimData$habitatGrid,
                              numGridRows = nrow(nimData$habitatGrid),
                              numGridCols = ncol(nimData$habitatGrid))
  }#i
  
  sex.init <- rbinom(n = nimConstants$n.individuals, 1, prob = 0.5)
  sex.init[!is.na(nimData$sex)] <- NA
  sex1.init <- sex.init + 1
  
  status.init <- rcat(n = nimConstants$n.individuals, prob = c(0.5,0.45,0.05))
  status.init[!is.na(nimData$status)] <- NA
  
  z.init <- rbinom(n = nimConstants$n.individuals, 1, prob = 0.1)
  z.init[!is.na(nimData$z)] <- NA
  
  nimInits <- list( "s" = s.init,
                    "z" = z.init,
                    "sex" = sex.init,
                    "status" = status.init,
                    "psi" = 0.5,
                    "rho" = 0.5,
                    "theta" = cbind(c(0.5,0.45,0.05),
                                    c(0.5,0.3,0.2)),
                    "betaDet" = rep(0,nimConstants$n.detCovs),
                    "betaHab" = rep(0,nimConstants$n.habCovs),
                    "p0" = cbind(c(0.1,0.1,0.05),
                                 c(0.1,0.1,0.05)),
                    "sigma" = cbind(c(1,1,2),
                                    c(1,1,2))
                    
                    save( modelCode,
                          nimData,
                          nimConstants,
                          nimInits,
                          nimParams,
                          file = file.path(thisDir, "input",
                                           paste0(modelName, "_", c, ".RData")))
} 


## ------   4. FIT MODEL -----
##---- Create the nimble model object
nimModel <- nimbleModel( code = modelCode,
                         constants = nimConstants,
                         inits = nimInits,
                         data = nimData,
                         check = FALSE,
                         calculate = FALSE) 
#nimModel$initializeInfo()
#nimModel$simulate()
nimModel$calculate()

CsimModel <- compileNimble(nimModel)
CsimModel$calculate()

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
                  bite.number = 20,
                  path = file.path(thisDir, paste0("output/chain",c))) 
  ))
}#c

##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 0,
                               param.omit = c("s","z"))

##---- Save MCMC output & plot
save( nimOutput,
      file = file.path(thisDir, paste0(modelName, "_mcmc.RData")))

pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput)
graphics.off()

plot(nimOutput)


## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   1. CALCULATE AC-BASED DENSITY ------
##-- Load & process MCMC samples
nimOutput_full <- collectMCMCbites( path = file.path(thisDir, "output"),
                                    burnin = 4)
res <- ProcessCodaOutput(nimOutput_full)

## Save processed MCMC samples
save(nimOutput_full, res,
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))

## LOAD AC-based DENSITY FILES
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))

## Extract the habitat raster
habitat.r <- habitat$raster

## Set cells outside the habitat to NA
habitat.r[habitat.r == 0] <- NA

## Create a matrix of raster cellsIDs
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = dim(habitat.r)[1],
                      ncol = dim(habitat.r)[2],
                      byrow = TRUE)

## Create a matrix of admin units binary indicators
## (rows == regions ; columns == habitat raster cells)
hab.rgmx <- rbind(habitat$raster[] == 1)
hab.rgmx[is.na(hab.rgmx)] <- 0
row.names(hab.rgmx) <- "habitat"

## Calculate density
WA_Density <- GetDensity_PD(
  sx = res$sims.list$s[ , ,1],
  sy = res$sims.list$s[ , ,2],
  z = res$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = hab.rgmx)


## Create a matrix of italy 
## (rows == regions ; columns == habitat raster cells)
it.rgmx <- rbind(habitat$Italia[] == 1)
it.rgmx[is.na(it.rgmx)] <- 0
row.names(it.rgmx) <- "Italia"

## Calculate density
WA_Italy <- GetDensity_PD(
  sx = res$sims.list$s[ , ,1],
  sy = res$sims.list$s[ , ,2],
  z = res$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = it.rgmx)


## Create a matrix of admin units binary indicators
## (rows == regions ; columns == habitat raster cells)
comp.rgmx <- rbind(habitat$comparison[] == 1)
comp.rgmx[is.na(comp.rgmx)] <- 0
row.names(comp.rgmx) <- "comparison"

## Calculate density
WA_Comp <- GetDensity_PD(
  sx = res$sims.list$s[ , ,1],
  sy = res$sims.list$s[ , ,2],
  z = res$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = T,
  regionID = comp.rgmx)




## ------   2. DENSITY MAPS ------
# meanDensity.R <- habitat.r
# meanDensity.R[] <- WA_Density$MeanCell
# meanDensity.R[is.na(habitat.r[])] <- NA
# writeRaster( x = meanDensity.R, overwrite = TRUE,
#              filename = file.path(thisDir, paste0(modelName, "_density.tif")))

pdf(file = file.path(thisDir, paste0(modelName,"_Density.pdf")),
    width = 20, height = 15)

##-- Set color scale
maxDens <- max(WA_Density$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)
par(mfrow = c(1,1), mar = c(4,4,0,0))

##-- Plot overall density raster
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


pdf(file = file.path(thisDir, paste0(modelName,"_IT_Density.pdf")),
    width = 15, height = 15)
##-- Plot Italian density raster
meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_Italy$MeanCell
meanDensity.R[is.na(habitat$Italia[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot( st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Italy$summary["Total",1],1),
                     " [", round(WA_Italy$summary["Total",4],1), " ; ",
                     round(WA_Italy$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)

graphics.off()
##-- Plot density raster for comparison between models
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





## ------  REALIZED vs. PREDICTED DENSITY  
##-- Calculate relative density
relativeDens.r <- meanDensity.R
relativeDens.r[] <- meanDensity.R[]/sum(meanDensity.R[], na.rm = T)

##-- Calculate relative point-process intensity
intensity.r <- habitat.r
intensity.r[intensity.r[] > 0] <- exp(res$mean$betaHab[1]*nimData$hab.covs[ ,1]+
                                        res$mean$betaHab[2]*nimData$hab.covs[ ,2]+
                                        res$mean$betaHab[3]*nimData$hab.covs[ ,3]+
                                        res$mean$betaHab[4]*nimData$hab.covs[ ,4]+
                                        res$mean$betaHab[5]*nimData$hab.covs[ ,5])
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

plot( habitat$polygon, border = grey(0.3))
mtext( text = "Predicted density", side = 3, line = -4, font = 2)
plot( relativeInt.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T)

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

centroids <- st_drop_geometry(ngs) %>%
  group_by(Genotype.ID) %>%
  summarise(X = mean(Coordinate.WGS84.UTM.East),               ## Get total length searched in each detector grid cell
            Y =  mean(Coordinate.WGS84.UTM.North)) 
coordinates(centroids) <- cbind.data.frame(centroids$X,
                                           centroids$Y)
proj4string(centroids) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
centroids <- st_as_sf(centroids)

plot( habitat$polygon, border = grey(0.3))
plot( relativeDens.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot(  habitat$polygon, add = T)
plot(st_geometry(countries), add = T, lwd=2)
plot(centroids,
     add = T, pch = 21, cex = 1.2,
     bg = adjustcolor("white",0.5),
     col = "white")
mtext( text = "Wolf density with\ndetected individuals' centroid",
       side = 3, font = 2, cex = 2)


graphics.off()





## -----------------------------------------------------------------------------
## ------   3. DETECTION ------
pdf(file = file.path(thisDir, paste0(myVars$modelName,"_Detection.pdf")),
    width = 20, height = 15)

par(mfrow = c(1,2))

## ------     3.1. DETECTION MAP ------
##-- DETECTABILITY MAP (p0) 
Detectors.list$main.detector.sp$p0 <- ilogit(logit(res$mean$p0[nimData$trapIntIndex]) +
                                               res$mean$betaTrap * nimData$DetsCovs)

pZero.r <- raster( extent(Detectors.list$main.detector.sp),
                   res = myVars$DETECTORS$detResolution)

pZero.r <- rasterize( x = Detectors.list$main.detector.sp,
                      y = pZero.r,
                      field = "p0")

##-- Set color scale
maxp0 <- round(max(Detectors.list$main.detector.sp$p0),3)
cuts <- seq(0, maxp0, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)

##-- plot p0 map
plot( Cleanpoly, col = "gray80", border = grey(0.3))
plot( pZero.r, add = T, box = F, bty = "n", 
      breaks = cuts, col = col,
      legend.width = 1,  
      axis.args = list( at = round(seq(0,maxp0,length.out = 5), digits = 5),
                        labels = round(seq(0,maxp0,length.out = 5), digits = 5), 
                        cex.axis = 0.6),
      legend.args = list( text = 'p0', line = -2,
                          side = 4, font = 2, cex = 0.8))



## ------     3.2. DETECTION EFFECT PLOT ------
par(mar = c(20,10,15,5))
detCovs <- seq(min(nimData$DetsCovs), max(nimData$DetsCovs), length.out = 100)
p0.posterior1 <- do.call( rbind,
                          lapply(1:length(res$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(res$sims.list$p0[x,1]) +
                                            res$sims.list$betaTrap[x] * detCovs)
                                 }))
p0.posterior2 <- do.call( rbind,
                          lapply(1:length(res$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(res$sims.list$p0[x,2]) +
                                            res$sims.list$betaTrap[x] * detCovs)
                                 }))
p0.posterior3 <- do.call( rbind,
                          lapply(1:length(res$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(res$sims.list$p0[x,3]) +
                                            res$sims.list$betaTrap[x] * detCovs)
                                 }))

maxP0 <- round(max(c(p0.posterior1,p0.posterior2,p0.posterior3)),3)    
p0cols <- hcl.colors(3)
plot( x = detCovs, y = colMeans(p0.posterior1),
      type = "n", ylim = c(0, maxP0), xlim = range(detCovs),
      ylab = "p0", xlab = "search tracks length (km)", axes = FALSE)

minTracks <- min(trackLength[])
maxTracks <- max(trackLength[])
xLabels <- round(seq(minTracks, maxTracks, length.out = 10),2)
axis(1, at = round(seq(min(detCovs), max(detCovs), length.out = 10),3),
     labels = xLabels, 
     tck = 0.01, las = 1, hadj = 0.5)
axis(2, at = seq(0,maxP0,length.out = 6),
     labels = seq(0,maxP0,length.out = 6),
     tck = 0.01, las = 1)

for(i in seq(1,80000,10)){
  points( x = detCovs, y = p0.posterior1[i, ], 
          type = "l", col = adjustcolor(p0cols[1], alpha.f = 0.01))
  points( x = detCovs, y = p0.posterior2[i, ], 
          type = "l", col = adjustcolor(p0cols[2], alpha.f = 0.01))
  points( x = detCovs, y = p0.posterior3[i, ], 
          type = "l", col = adjustcolor(p0cols[3], alpha.f = 0.01))
}

points( x = detCovs, y = colMeans(p0.posterior1), lwd = 2, 
        type = "l", col = p0cols[1])
points( x = detCovs, y = colMeans(p0.posterior2), lwd = 2, 
        type = "l", col = p0cols[2])
points( x = detCovs, y = colMeans(p0.posterior3), lwd = 2, 
        type = "l", col = p0cols[3])

dev.off()



## ------   4. PARAMETER TABLES ------
params.simple <- names(res$mean)[2:4]
params.means <- do.call(c, lapply(params.simple, function(x)res$mean[[x]][1]))
params.sd <- do.call(c, lapply(params.simple, function(x)res$sd[[x]][1]))
params.q2.5 <- do.call(c, lapply(params.simple, function(x)res$q2.5[[x]][1]))
params.q97.5 <- do.call(c, lapply(params.simple, function(x)res$q97.5[[x]][1]))
params.Rhat <- do.call(c, lapply(params.simple, function(x)res$Rhat[[x]][1]))
params.n.eff <- do.call(c, lapply(params.simple, function(x)res$n.eff[[x]][1]))

params.summary <- cbind.data.frame( params.simple, params.means, params.sd,
                                    params.q2.5, params.q97.5, params.Rhat,
                                    params.n.eff)

names(params.summary) <- c("params", "mean", "sd", "2.5%CI", "97.5%CI", "Rhat", "n.eff")

## WRITE TABLE
write.csv( x = params.summary,
           file = file.path(myVars$WD, myVars$modelName, "TABLES",
                            paste( myVars$modelName, "_params.csv", sep = "")))

print( xtable(params.summary, type = "latex"),
       floating = FALSE,# scalebox=.8,
       add.to.row = list(list(seq(1, nrow(params.summary), by = 2)),"\\rowcolor[gray]{.95} "),
       file = file.path(myVars$WD, myVars$modelName,"TABLES",
                        paste( myVars$modelName, "_params.tex", sep = "")))


## -----------------------------------------------------------------------------
{## ------ IV. EXPORT RESULTS REPORT AS .pdf ------
  ##-- Load PROCESSED MODEL OUTPUTS
  load(file.path(myVars$WD, myVars$outFolder, "results_F.RData"))
  myResults_F <- myResults
  nimOutput.simple_F <- nimOutput.simple
  rm("myResults", "nimOutput.simple")
  
  load(file.path(myVars$WD, myVars$outFolder, "results_M.RData"))
  myResults_M <- myResults
  nimOutput.simple_M <- nimOutput.simple
  rm("myResults", "nimOutput.simple")
  
  
  ## LOAD PRE_CALCULATED DENSITIES 
  load(file.path(myVars$WD, myVars$outFolder, "Density_admin.RData" ))
  load(file.path(myVars$WD, myVars$outFolder, "Density_country.RData" ))
  
  ## LOAD ONE INPUT FILE FOR EACH SEX
  inFiles <- list.files(file.path(myVars$WD, myVars$modelNameF, "NimbleInFiles"))
  load(file.path(myVars$WD, myVars$modelNameF, "NimbleInFiles", inFiles[1]))
  inits_F <- nimInits
  
  inFiles <- list.files(file.path(myVars$WD, myVars$modelNameM, "NimbleInFiles"))
  load(file.path(myVars$WD, myVars$modelNameM, "NimbleInFiles", inFiles[1]))
  inits_M <- nimInits
  
  ## LOAD HABITAT
  load(file.path(dir.dropboxALL, "BavRedDeer/HABITAT.RData"))
  
  ## LOAD DETECTORS
  load(file = file.path(dir.dropboxALL,"BavRedDeer/DETECTORS.RData"))
  
  ## LOAD NGS DATA AND SHAPE FILES
  load(file.path(myVars$WD, myVars$outFolder, "processedObjects.RData" ))
  
  
  
  ################### START PDF ################# 
  pdf( file = file.path(myVars$WD, myVars$outFolder, "ResultsCompilation2.pdf"),
       width = 20, height = 12)
  ###############################################
  ## ------   1. BEFORE/AFTER PLOTS ------
  # pdf( file = file.path(myVars$WD, myVars$outFolder, "PLOTS",
  #                       paste(myVars$outFolder, "_data_prep.pdf", sep = "")),
  #      width = 20, height = 12)
  
  
  ## ------     1.1. DATA PROCESSING ------
  par(mfrow = c(1,3))
  ##-- RAW DATA
  ##-- Study area polygon with admin units
  plot(Cleanpoly, col = adjustcolor(hcl.colors(n = 3), alpha.f = 0.5), border = F)
  plot(Cleanpoly, col = adjustcolor("black", alpha.f = 0.3), border = F, add = T)
  mtext(text = "Search tracks\nand samples", line = -8, side = 3, font = 2, cex = 1.5)
  ##-- GPS search tracks 
  plot(Tracks, col = "firebrick3", add = T, lwd = 0.5)
  ##-- NGS samples per sex
  sexCol <- ifelse(Data$Sex == "F", "black", "gray80")
  sexCol[is.na(sexCol)] <- "gray40"
  plot(myData.sp$myData.sp, add = T, pch = 19, col = sexCol, cex = 0.5)
  legend( x = 372000, y = 5410000, bty = "n",
          legend = c("Females", "Males"), title = "NGS samples",
          pch = 19, col = c("black", "gray80"))
  
  
  ##-- PROCESSED DATA
  ##-- NUMBER OF NGS SAMPLES
  plot(Cleanpoly, col = adjustcolor(hcl.colors(n = 3), alpha.f = 0.5), border = F)
  mtext( text = "Number of samples per detector", line = -8, side = 3, font = 2, cex = 1.5)
  ##-- Habitat raster
  temp <- Habitat.list$habitat.r
  temp[temp[] == 0] <- NA
  plot(temp, add = T, col = adjustcolor("black", alpha.f = 0.3), legend = F)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  ##-- Main detectors
  plot(Detectors.list$main.detector.sp, add = T, cex = 0.3, col = "gray40")
  ##-- Number of detections per detector
  numDets <- apply(y, 2, sum)
  colPal <- hcl.colors(n = max(numDets)+1)
  points(Detectors.list$main.detector.sp,
         cex = numDets/2,
         pch = 21, bg = colPal[numDets+1])
  
  ##-- NUMBER OF IDs DETECTED
  plot(Cleanpoly, col = adjustcolor(hcl.colors(n = 3), alpha.f = 0.5), border = F)
  mtext( text = "Number of individuals per detector", line = -8, side = 3, font = 2, cex = 1.5)
  ##-- Habitat raster
  temp <- Habitat.list$habitat.r
  temp[temp[] == 0] <- NA
  plot(temp, add = T, col = adjustcolor("black", alpha.f = 0.3), legend = F)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  ##-- Main detectors
  plot(Detectors.list$main.detector.sp, add = T, cex = 0.3, col = "gray40")
  ##-- Number of individuals per detector
  numIdDets <- apply(y, 2, function(x)sum(x>0))
  points(Detectors.list$main.detector.sp,
         cex = numIdDets/2,
         pch = 21, bg = colPal[numIdDets+1])
  
  mtext( text = "Data Processing", side = 3, line = -5, outer = TRUE, font = 2, cex = 2)
  
  levels <- unique(as.numeric(c(names(table(numIdDets)), names(table(numDets)))))
  nameNumDets <- levels[order(levels)][-1]
  legend( x = 372000, y = 5410000, bty = "n",
          legend = rev(nameNumDets), title = "# of obs/ids per detector",
          y.intersp = 1.2, x.intersp = 2, pch = 21,
          pt.bg = rev(colPal[nameNumDets+1]),
          pt.cex = rev((nameNumDets+1)/2))
  
  
  
  ## ------     1.2. SPATIAL COVARIATES ------
  ##-- DISTURBANCES 
  par(mfrow = c(2,3))
  ##-- Raw data
  cuts <- 0:4  
  cols <- hcl.colors(4) #rev(terrain.colors(4))
  ##-- All disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "All disturbances (raw)", side = 1, font = 2)
  temp <- Dist_agents
  temp[!Dist_year[] %in% c(1995:2020)] <- NA
  plot( temp, add = T,
        breaks = cuts, col = cols,
        axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  ##-- Managed Disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Managed disturbances (raw)", side = 1, font = 2)
  temp <- Dist_agents
  temp[!(Dist_year[] %in% c(1995:2020) & Dist_agents[] %in% c(3,4))] <- NA
  plot(temp , add = T, 
       breaks = cuts, col = cols,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  ##-- Unmanaged Disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Unmanaged disturbances (raw)", side = 1, font = 2)
  temp <- Dist_agents
  temp[!(Dist_year[] %in% c(1995:2020) & Dist_agents[] %in% c(1,2))] <- NA
  plot(temp , add = T,   
       breaks = cuts, col = cols,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  
  ##-- Processed data
  ##-- All disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "All disturbances (processed)", side = 1, font = 2)
  plot( disturb.r, add = T,
        axes = F, box = F, bty = "n", legend = F)
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), border = "gray40", add =T)
  plot( Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  ##-- Managed Disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Managed disturbances (processed)", side = 1, font = 2)
  plot( managed.r, add = T,
        axes = F, box = F, bty = "n", legend = F)
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), border = "gray40", add =T)
  plot( Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  ##-- Unmanaged Disturbances
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Unmanaged disturbances (processed)", side = 1, font = 2)
  plot( unmanaged.r, add = T,
        axes = F, box = F, bty = "n", legend = F)
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), border = "gray40", add =T)
  plot( Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  mtext( text = "Spatial covariates", side = 3, line = -3, cex = 2, outer = TRUE, font = 2)
  
  
  par(mfrow = c(2,2), mar = c(2,2,2,2))
  ##-- Raw data
  ##-- ELEVATION
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Elevation (raw)", side = 1, font = 2)
  plot( Elevation.r, add = T, 
        axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add =T)
  
  ##-- Search tracks
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3))
  mtext( text = "Search tracks (raw)", side = 1, font = 2)
  plot(Tracks, col = "firebrick3", add = T, lwd = 0.5)
  plot(Cleanpoly, add = T)
  
  
  ##-- Processed data
  ##-- ELEVATION
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Elevation (processed)", side = 1, font = 2)
  plot( DEM.r, add = T, 
        axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), border = "gray40", add =T)
  plot( Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  mtext( text = "Spatial covariates", side = 3, line = -3, outer = TRUE, font = 2, cex = 2)
  
  ##-- Search-track length
  plot( Cleanpoly, border = grey(0.95))
  Detectors.list$main.detector.sp$trackLength <- nimData$DetsCovs
  r <- rasterFromXYZ(as.data.frame(Detectors.list$main.detector.sp)[ ,c("x","y","trackLength")])
  plot( r, add = T, axes = F, box = F, bty = "n", legend = F)
  
  mtext( text = "Search tracks length (processed)", side = 1, font = 2)
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), border = "gray40", add =T)
  plot( Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  
  ##-- REGIONS 
  par(mfrow = c(2,3), mar = c(2,2,2,2))
  ##-- All regions
  cuts <- 0:5  
  cols <- hcl.colors(5) # rev(terrain.colors(4))
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "Regions", side = 1, font = 2)
  plot( REGIONS.r, add = T,
        breaks = cuts, col = cols,
        axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  
  ##-- BFNP
  cuts <- seq(0, 1, length.out = 11)   
  colFunc <- colorRampPalette(c("#FFFFFF", cols[1]))
  col <- colFunc(10)
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "BFNP", side = 1, font = 2)
  plot(REGION_BNF.r, add = T,
       breaks = cuts, col = col,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  ##-- SUM
  cuts <- seq(0, 1, length.out = 11)   
  colFunc <- colorRampPalette(c("#FFFFFF", cols[3]))
  col <- colFunc(10)
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "SNP", side = 1, font = 2)
  plot(REGION_SUM.r, add = T,   
       breaks = cuts, col = col,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  ##-- NR
  cuts <- seq(0, 1, length.out = 11)   
  colFunc <- colorRampPalette(c("#FFFFFF", cols[2]))
  col <- colFunc(10)
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "NR", side = 1, font = 2)
  plot(REGION_NEU.r, add = T,   
       breaks = cuts, col = col,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  ##-- BFNP (protected)
  cuts <- seq(0, 1, length.out = 11)   
  colFunc <- colorRampPalette(c("#FFFFFF", cols[4]))
  col <- colFunc(10)
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "BFNP (protected)", side = 1, font = 2)
  plot(REGION_BNFp.r, add = T,   
       breaks = cuts, col = col,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  
  ##-- SUM (protected)
  cuts <- seq(0, 1, length.out = 11)   
  colFunc <- colorRampPalette(c("#FFFFFF", cols[5]))
  col <- colFunc(10)
  plot( Cleanpoly, border = grey(0.95))
  mtext( text = "SNP (protected)", side = 1, font = 2)
  plot(REGION_SUMp.r, add = T,   
       breaks = cuts, col = col,
       axes = F, box = F, bty = "n")
  plot(Cleanpoly, col = adjustcolor("gray60", alpha.f = 0.3), add = T)
  plot(Habitat.list$buffered.habitat.poly, add = T, border = "red", lwd = 2)
  
  
  mtext( text = "Spatial covariates", side = 3, line = -3, cex = 2, outer = TRUE, font = 2)
  
  # dev.off()
  
  
  ## ------   2. DENSITY PLOTS ------
  # pdf(file = file.path( myVars$WD, myVars$outFolder, "PLOTS",
  #                       paste(myVars$outFolder, "_Density.pdf", sep = "")),
  #     width = 20, height = 15)
  
  ## ------     2.1. DENSITY MAPS ------
  ##-- Set color scale
  maxDens <- max(c(Density$MeanCell, SpaceUse$MeanCell))
  cuts <- seq(0, maxDens, length.out = 100)   
  colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
  #colFunc <- colorRampPalette(c("white","red"))
  col <- colFunc(100)
  
  
  ## ------       2.1.a. OVERALL DENSITY MAP ------
  habitat.r <- Habitat.list$habitat.r
  habitat.r[habitat.r == 0] <- NA
  admin.r <- rasterize(Cleanpoly, habitat.r)
  
  meanDensity.R <- habitat.r
  meanDensity.R[ ] <- Density$MeanCell
  meanDensity.R[is.na(admin.r[])|is.na(habitat.r[])] <- NA
  
  par(mfrow = c(1,1), mar = c(4,4,0,0))
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  mtext( text = paste( "N = ", round(Density$summary["Total",1],1),
                       " [", round(Density$summary["Total",4],1), " ; ",
                       round(Density$summary["Total",5],1), "]", sep = ""),
         side = 1, font = 2, cex = 1.5)
  mtext( text = "Total density", side = 3, font = 2, cex = 2, line = -2)
  
  ##-- Add region-specific pop size
  center.region <- rbind( c(382330,#mean(extent(Cleanpoly_BFNP)[1:2]),
                            5423094),# mean(extent(Cleanpoly_BFNP)[3:4])),
                          c(405084,#mean(extent(Cleanpoly_NR)[1:2]),
                            5410600),#mean(extent(Cleanpoly_NR)[3:4])),
                          c(398095,#mean(extent(Cleanpoly_SNP)[1:2]),
                            5428000))#mean(extent(Cleanpoly_SNP)[3:4])))
  center.region <- SpatialPoints(coords = center.region)
  # text( x = center.region,
  #       labels = round(Density$summary[Cleanpoly$Unit,1],1),
  #       cex = 1.5, font = 2)
  
  ##-- Add density legend
  plot( meanDensity.R, breaks = cuts, col = col,
        legend.width = 1, legend.only = TRUE,
        axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
                          labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
                          cex.axis = 0.6),
        legend.args = list( text = 'Density (ids.km2)', line = -2,
                            side = 4, font = 2, cex = 0.8))
  
  
  
  
  ## ------       2.1.c. OVERALL SPACE-USE MAP ------
  habitat.r <- Habitat.list$habitat.r
  habitat.r[habitat.r == 0] <- NA
  admin.r <- rasterize(Cleanpoly, habitat.r)
  
  meanDensity.R <- habitat.r
  meanDensity.R[!is.na(meanDensity.R[])] <- SpaceUse$MeanCell
  meanDensity.R[is.na(admin.r[])|is.na(habitat.r[])] <- NA
  
  par(mfrow = c(1,1), mar = c(4,4,0,0))
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  mtext( text = paste( "N = ", round(SpaceUse$summary["Total",1],1),
                       " [", round(SpaceUse$summary["Total",4],1), " ; ",
                       round(SpaceUse$summary["Total",5],1), "]", sep = ""),
         side = 1, font = 2, cex = 1.5)
  mtext( text = "Total utilization distribution", side = 3, font = 2, cex = 2, line = -2)
  
  ##-- Add region-specific pop size
  center.region <- rbind( c(382330,#mean(extent(Cleanpoly_BFNP)[1:2]),
                            5423094),# mean(extent(Cleanpoly_BFNP)[3:4])),
                          c(405084,#mean(extent(Cleanpoly_NR)[1:2]),
                            5410600),#mean(extent(Cleanpoly_NR)[3:4])),
                          c(398095,#mean(extent(Cleanpoly_SNP)[1:2]),
                            5428000))#mean(extent(Cleanpoly_SNP)[3:4])))
  center.region <- SpatialPoints(coords = center.region)
  # text( x = center.region,
  #       labels = round(SpaceUse$summary[Cleanpoly$Unit,1],1),
  #       cex = 1.5, font = 2)
  
  ##-- Add density legend
  plot( meanDensity.R, breaks = cuts, col = col,
        legend.width = 1, legend.only = TRUE,
        axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
                          labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
                          cex.axis = 0.6),
        legend.args = list( text = 'Density (ids.km2)', line = -2,
                            side = 4, font = 2, cex = 0.8))
  
  
  
  
  
  
  ## ------       2.1.b. SEX-SPECIFIC DENSITY MAPS ------
  par(mfrow = c(1,2), mar = c(4,4,0,0))
  ##-- FEMALES
  meanDensity_F.R <- habitat.r
  meanDensity_F.R[ ] <- Density_F$MeanCell
  meanDensity_F.R[is.na(admin.r[])|is.na(habitat.r[])] <- NA
  
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  mtext( text = "Female density", side = 3, font = 2, cex = 2, line = -2)
  plot( meanDensity_F.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  mtext( text = paste( "N = ", round(Density_F$summary["Total",1],1),
                       " [", round(Density_F$summary["Total",4],1), " ; ",
                       round(Density_F$summary["Total",5],1), "]", sep = ""),
         side = 1, font = 2, cex = 1.5)
  ##-- Add region-specific pop size
  # text( x = center.region,
  #       labels = round(Density_F$summary[Cleanpoly$Unit,1],1),
  #       cex = 1.2, font = 2)
  
  ##-- MALES
  par(mar = c(4,0,0,4))
  meanDensity_M.R <- habitat.r
  meanDensity_M.R[ ] <- Density_M$MeanCell
  meanDensity_M.R[is.na(admin.r[])|is.na(habitat.r[])] <- NA
  
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  mtext( text = "Male density", side = 3, font = 2, cex = 2, line = -2)
  plot( meanDensity_M.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  mtext( text = paste( "N = ", round(Density_M$summary["Total",1],1),
                       " [", round(Density_M$summary["Total",4],1), " ; ",
                       round(Density_M$summary["Total",5],1), "]", sep = ""),
         side = 1, font = 2, cex = 1.5)
  
  ##-- Add region-specific pop size
  # text( x = center.region,
  #       labels = round(Density_M$summary[Cleanpoly$Unit,1],1),
  #       cex = 1.2, font = 2)
  
  ## Add density legend
  plot( meanDensity.R, breaks = cuts, col = col,
        legend.width = 1, legend.only = TRUE,
        axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
                          labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
                          cex.axis = 0.6),
        legend.args = list( text = 'Density (ids.km2)', line = -2,
                            side = 4, font = 2, cex = 0.8))
  
  
  ## ------       2.1.d. DETECTION CENTROIDS PLOT ----- 
  par(mfrow = c(1,2), mar = c(4,4,0,0))
  
  ##-- FEMALES
  individual.centroids <- inits_F$sxy[is.na(inits_F$z), ]
  colnames(individual.centroids) <- c("x","y")
  hab.xy <- coordinates(Habitat.list$habitat.sp)
  colnames(hab.xy) <- c("x","y")
  individual.centroids <- scaleCoordsToHabitatGrid(coordsData = individual.centroids,
                                                   scaleToGrid = FALSE,
                                                   coordsHabitatGridCenter = hab.xy)
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  mtext( text = "Female density with\nindividual detection centroids",
         side = 3, font = 2, cex = 2, line = -4)
  plot( meanDensity_F.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  points( individual.centroids$coordsDataScaled, cex = 1.5,
          pch = 19, col = adjustcolor("black", alpha.f = 0.2))
  
  
  ##-- MALES
  individual.centroids <- inits_M$sxy[is.na(inits_M$z), ]
  colnames(individual.centroids) <- c("x","y")
  hab.xy <- coordinates(Habitat.list$habitat.sp)
  colnames(hab.xy) <- c("x","y")
  individual.centroids <- scaleCoordsToHabitatGrid(coordsData = individual.centroids,
                                                   scaleToGrid = FALSE,
                                                   coordsHabitatGridCenter = hab.xy)
  par(mar = c(4,0,0,4))
  plot( Cleanpoly, col = "gray80", border = grey(0.3))
  mtext( text = "Male density with\nindividual detection centroids",
         side = 3, font = 2, cex = 2, line = -4)
  plot( meanDensity_M.R, add = T,
        breaks = cuts, col = col, colNA = NA,
        axes = F, box = F, bty = "n", legend = F)
  plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
  points( individual.centroids$coordsDataScaled, cex = 1.5,
          pch = 19, col = adjustcolor("black", alpha.f = 0.2))
  
  # dev.off()
  
  
  ## ------     2.2. DENSITY EFFECT PLOT ------
  ##-- Check the range of covariate values
  range(nimData$DEM)
  range(nimData$R)
  
  ##-- Create a data.frame of covariate values
  data <- cbind.data.frame( elev = as.vector(as.character(round(nimData$DEM, digits = 1))), 
                            dist = as.vector(as.character(round(nimData$R, digits = 2))))
  data <- unique(data)
  data[data == "-1"] <- -1.0
  data[data == "1"] <- 1.0
  data[data == "0"] <- 0.0
  data$test <- paste(data$elev, data$dist)
  
  ##-- PREDICTED DENSITY FOR MEAN "unmanaged disturbances" :
  ##-- Generate covariate values to predict sigma and HR for
  ##-- Weird way of generating values to avoid rounding problems later and 
  ##-- generate sigma and HR values only for the covariate values encountered in real life ...
  ##-- It's ugly and there is probably a simple fix for that :-)
  elev <- (-22:24) * 0.1
  dist <- (0:90) * 0.01
  predDensity_F <- predDensity_M <- matrix(NA, nrow = length(elev), ncol = length(dist))
  for(e in 1:length(elev)){
    for(l in 1:length(dist)){
      if(paste(elev[e],dist[l]) %in% data$test){ 
        predDensity_F[e,l] <- exp(
          myResults_F$mean$betaR * dist[l] +
            myResults_F$mean$betaDEM * elev[e] +
            myResults_F$mean$betaDEMR * dist[l] * elev[e])
        predDensity_M[e,l] <- exp(
          myResults_M$mean$betaR * dist[l] +
            myResults_M$mean$betaDEM * elev[e] +
            myResults_M$mean$betaDEMR * dist[l] * elev[e]) 
      }
    }
  }
  
  ##-- 3D PLOTS
  par(mfrow= c(1,2), mar = c(2,2,0,0))
  ##-- plot the box in the background
  perspbox(x = elev, y = dist, z = predDensity_F,
           xlab = "Elevation", ylab = "All disturbances", zlab = "Female density",
           zlim = c(0,6),
           bty = "u",
           col.panel = "gray80", ## BACKGROUND COLOR
           col.axis = "black",   ## GRID COLOR
           col.grid = "white",   ## AXES COLOR
           ticktype = "detailed",
           theta = -45,           ## ANGLE FOR ORIENTATION LEFT/RIGHT
           phi = 20)             ## ANGLE FOR ORIENTATION UP/DOWN
  ##-- plot the histogram on top
  hist3D(x = elev, y = dist, z = predDensity_F, scale = TRUE,
         col = ramp.col(c("blue", "yellow", "firebrick")), #terrain.colors(10)
         border = "gray40",
         colkey = list( side = 1, plot = F),
         shade = 0.3,
         space = 0.0,
         add = T)
  
  
  ##-- plot the box in the background
  perspbox(x = elev, y = dist, z = predDensity_M,
           xlab = "Elevation", ylab = "All disturbances", zlab = "Male density",
           zlim = c(0,6),
           bty = "u",
           col.panel = "gray80", ## BACKGROUND COLOR
           col.axis = "black",   ## GRID COLOR
           col.grid = "white",   ## AXES COLOR
           ticktype = "detailed",
           theta = -45,           ## ANGLE FOR ORIENTATION LEFT/RIGHT
           phi = 20)             ## ANGLE FOR ORIENTATION UP/DOWN
  ##-- plot the histogram on top
  hist3D(x = elev, y = dist, z = predDensity_M, scale = TRUE,
         col = ramp.col(c("blue", "yellow", "firebrick")), #terrain.colors(10)
         border = "gray40",
         colkey = list( side = 1, plot = F),
         shade = 0.3,
         space = 0.0,
         add = T)
  # dev.off()
  
  
  
  ## ------   3. DENSITY TABLES ------
  temp <- admin.r
  temp[is.na(habitat.r[])] <- NA 
  numCells_SNP <- sum(temp[] %in% 3)
  numCells_BFNP <- sum(temp[] %in% 1)
  numCells_NR <- sum(temp[] %in% 2)
  numCells_GER <- numCells_BFNP + numCells_NR
  numCells_Total <- numCells_BFNP + numCells_NR + numCells_SNP
  numCells <- c(numCells_SNP,numCells_BFNP, numCells_NR, numCells_GER, numCells_Total)
  
  ## FEMALE ABUNDANCES 
  mean_N_Females <- round(c(Density_F$summary[c(3,1,2),1],
                            Density_country_F$summary[c(1,3),1]),1)
  lower_N_Females <- c(Density_F$summary[c(3,1,2),4],
                       Density_country_F$summary[c(1,3),4])
  upper_N_Females <- c(Density_F$summary[c(3,1,2),5],
                       Density_country_F$summary[c(1,3),5])
  N_Females <- paste( mean_N_Females, " [", lower_N_Females,
                      ";",  upper_N_Females,"]", sep = "")
  
  ## FEMALE DENSITIES (ids/1)
  D_Females <- paste( round( mean_N_Females/numCells, 2), " [",
                      round( lower_N_Females/numCells, 2), ";",
                      round( upper_N_Females/numCells, 2), "]", sep = "")
  
  
  ## MALE ABUNDANCES 
  mean_N_Males <- round(c(Density_M$summary[c(3,1,2),1],
                          Density_country_M$summary[c(1,3),1]),1)
  lower_N_Males <- c(Density_M$summary[c(3,1,2),4],
                     Density_country_M$summary[c(1,3),4])
  upper_N_Males <- c(Density_M$summary[c(3,1,2),5],
                     Density_country_M$summary[c(1,3),5])
  N_Males <- paste( mean_N_Males, " [", lower_N_Males,
                    ";",  upper_N_Males,"]", sep = "")
  
  ## MALE DENSITIES (ids/1)
  mean_N_Males/numCells
  D_Males <- paste( round( mean_N_Males/numCells, 2), " [",
                    round( lower_N_Males/numCells, 2), ";",
                    round( upper_N_Males/numCells, 2), "]", sep = "")
  
  
  
  ## TOTAL ABUNDANCES 
  mean_N_Total <- round(c(Density$summary[c(3,1,2),1],
                          Density_country$summary[c(1,3),1]),1)
  lower_N_Total <- c(Density$summary[c(3,1,2),4],
                     Density_country$summary[c(1,3),4])
  upper_N_Total <- c(Density$summary[c(3,1,2),5],
                     Density_country$summary[c(1,3),5])
  N_Total <- paste( mean_N_Total, " [", lower_N_Total,
                    ";",  upper_N_Total,"]", sep = "")
  
  mean_NSU_Total <- round(c(SpaceUse$summary[c(3,1,2),1],
                            sum(SpaceUse$summary[c(1,2),1]),
                            SpaceUse$summary[c(4),1]),1)
  lower_NSU_Total <- round(c(SpaceUse$summary[c(3,1,2),4],
                             sum(SpaceUse$summary[c(1,2),4]),
                             SpaceUse$summary[c(4),4]),1)
  upper_NSU_Total <- round(c(SpaceUse$summary[c(3,1,2),5],
                             sum(SpaceUse$summary[c(1,2),5]),
                             SpaceUse$summary[c(4),5]),1)
  NSU_Total <- paste( mean_NSU_Total, " [", lower_NSU_Total,
                      ";",  upper_NSU_Total,"]", sep = "")
  ## TOTAL DENSITIES (ids/1)
  SU_Total <- paste( round( mean_NSU_Total/numCells, 2), " [",
                     round( lower_NSU_Total/numCells, 2), ";",
                     round( upper_NSU_Total/numCells, 2), "]", sep = "")
  
  
  N_Table <- cbind.data.frame(N_Females, N_Males, N_Total, NSU_Total)
  names(N_Table) <- c("Females", "Males", "Total", "Total (UD)")
  row.names(N_Table) <- c("SNP","BFNP","NR","GER","Total")
  
  D_Table <- cbind.data.frame(D_Females, D_Males, D_Total, SU_Total)
  names(D_Table) <- c("Females", "Males", "Total", "Total (UD)")
  row.names(D_Table) <- c("SNP","BFNP","NR","GER","Total")
  
  
  ##-- Export as .pdf
  par(mfrow = c(1,1))
  plot.new()
  grid.table(N_Table)
  mtext(text = "Abundance estimates\nby region and sex",
        side = 3, outer = T, line = -20, cex = 2, font = 2)
  
  plot.new()
  grid.table(D_Table)
  mtext(text = "Mean density estimates\nby region and sex (ids/100ha) ",
        side = 3, outer = T, line = -20, cex = 2, font = 2)
  
  # ##-- Export as .csv
  # write.csv( x = N_Table,
  #            file = file.path(myVars$WD, myVars$outFolder, "TABLES",
  #                             paste( myVars$outFolder, "_N.csv", sep = "")))
  # ##-- Export as .tex
  # print( xtable(N_Table, type = "latex"),
  #        floating = FALSE, # scalebox=.8,
  #        add.to.row = list(list(seq(1, nrow(N_Table), by = 2)),"\\rowcolor[gray]{.95} "),
  #        file = file.path(myVars$WD, myVars$outFolder, "TABLES",
  #                         paste(myVars$outFolder, "_N.tex", sep = "")))
  
  
  ## ------   4. DETECTION ------
  # pdf(file = file.path(myVars$WD, myVars$outFolder, "PLOTS",
  #                      paste(myVars$outFolder,"_Detection.pdf", sep ="")),
  #     width = 20, height = 15)
  
  par(mfrow = c(1,2))
  
  
  # ## ------     4.1. DETECTION MAPS ------
  # ## ------       4.1.a. FEMALES ------
  # ##-- DETECTABILITY MAP (p0) 
  # Detectors.list$main.detector.sp$p0 <- ilogit(logit(myResults_F$mean$p0[1]) +
  #                                                myResults_F$mean$betaTrap * nimData$DetsCovs)
  # 
  # pZero.r <- raster( extent(Detectors.list$main.detector.sp),
  #                    res = myVars$DETECTORS$detResolution)
  # 
  # pZero.r <- rasterize( x = Detectors.list$main.detector.sp,
  #                       y = pZero.r,
  #                       field = "p0")
  # 
  # ##-- Set color scale
  # maxp0 <- round(max(Detectors.list$main.detector.sp$p0), 3)
  # cuts <- seq(0, maxp0, length.out = 100)   
  # colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
  # col <- colFunc(100)
  # 
  # ##-- plot p0 map
  # par(mar= c(4,4,0,0))
  # plot( Cleanpoly, col = "gray80", border = grey(0.3))
  # plot( pZero.r, add = T, box = F, bty = "n", 
  #       breaks = cuts, col = col,
  #       legend.width = 1,  
  #       axis.args = list( at = round(seq(0,maxp0,length.out = 5), digits = 5),
  #                         labels = round(seq(0,maxp0,length.out = 5), digits = 5), 
  #                         cex.axis = 0.6),
  #       legend.args = list( text = 'p0', line = -2,
  #                           side = 4, font = 2, cex = 0.8))
  # 
  # 
  # ## ------       4.1.b. MALES ------
  # ##-- DETECTABILITY MAP (p0) 
  # Detectors.list$main.detector.sp$p0 <- ilogit(logit(myResults_M$mean$p0[1]) +
  #                                                myResults_M$mean$betaTrap * nimData$DetsCovs)
  # 
  # pZero.r <- raster( extent(Detectors.list$main.detector.sp),
  #                    res = myVars$DETECTORS$detResolution)
  # 
  # pZero.r <- rasterize( x = Detectors.list$main.detector.sp,
  #                       y = pZero.r,
  #                       field = "p0")
  # 
  # ##-- Set color scale
  # maxp0 <- round(max(Detectors.list$main.detector.sp$p0),3)
  # cuts <- seq(0, maxp0, length.out = 100)   
  # colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
  # col <- colFunc(100)
  # 
  # ##-- plot p0 map
  # par(mar= c(4,0,0,4))
  # plot( Cleanpoly, col = "gray80", border = grey(0.3))
  # plot( pZero.r, add = T, box = F, bty = "n", 
  #       breaks = cuts, col = col,
  #       legend.width = 1,  
  #       axis.args = list( at = round(seq(0,maxp0,length.out = 5), digits = 5),
  #                         labels = round(seq(0,maxp0,length.out = 5), digits = 5), 
  #                         cex.axis = 0.6),
  #       legend.args = list( text = 'p0', line = -2,
  #                           side = 4, font = 2, cex = 0.8))
  
  
  
  ## ------     4.2. DETECTION EFFECT PLOT ------
  ## ------       4.2.a. FEMALES ------
  par(mar = c(20,10,15,5))
  detCovs <- seq(min(nimData$DetsCovs), max(nimData$DetsCovs), length.out = 100)
  p0.posterior1 <- do.call( rbind,
                            lapply(1:length(myResults_F$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_F$sims.list$p0[x,1]) +
                                              myResults_F$sims.list$betaTrap[x] * detCovs)
                                   }))
  p0.posterior2 <- do.call( rbind,
                            lapply(1:length(myResults_F$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_F$sims.list$p0[x,2]) +
                                              myResults_F$sims.list$betaTrap[x] * detCovs)
                                   }))
  p0.posterior3 <- do.call( rbind,
                            lapply(1:length(myResults_F$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_F$sims.list$p0[x,3]) +
                                              myResults_F$sims.list$betaTrap[x] * detCovs)
                                   }))
  
  maxP0 <- round(max(c(p0.posterior1,p0.posterior2,p0.posterior3)),3)    
  p0cols <- c("#4B0055", "#009B85", "#FDE333") #hcl.colors(3)
  plot( x = detCovs, y = colMeans(p0.posterior1),
        type = "n", ylim = c(0, maxP0), xlim = range(detCovs),
        ylab = "p0", xlab = "search tracks length (km)", axes = FALSE)
  
  minTracks <- min(trackLength[])
  maxTracks <- max(trackLength[])
  xLabels <- round(seq(minTracks, maxTracks, length.out = 10),2)
  axis(1, at = round(seq(min(detCovs), max(detCovs), length.out = 10),3),
       labels = xLabels, 
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxP0,length.out = 6),
       labels = seq(0,maxP0,length.out = 6),
       tck = 0.01, las = 1)
  
  for(i in seq(1,dim(p0.posterior1)[1],10)){
    points( x = detCovs, y = p0.posterior1[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[1], alpha.f = 0.01))
    points( x = detCovs, y = p0.posterior2[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[2], alpha.f = 0.01))
    points( x = detCovs, y = p0.posterior3[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[3], alpha.f = 0.01))
  }
  
  points( x = detCovs, y = colMeans(p0.posterior1), lwd = 2, 
          type = "l", col = p0cols[1])
  points( x = detCovs, y = colMeans(p0.posterior2), lwd = 2, 
          type = "l", col = p0cols[2])
  points( x = detCovs, y = colMeans(p0.posterior3), lwd = 2, 
          type = "l", col = p0cols[3])
  
  mtext( text = "Females", side = 3, line = -3, font = 2, cex = 1.5)
  
  mtext( text = "Effect of search track length\non baseline detection probability",
         side = 3, line = -6, outer = TRUE, font = 2, cex = 2)
  
  legend("topright", bty = "n",
         legend = c("BFNP","NR", "SNP"), 
         lwd = 2, col = p0cols)
  
  
  ## ------       4.2.b. MALES ------
  par(mar = c(20,10,15,5))
  detCovs <- seq(min(nimData$DetsCovs), max(nimData$DetsCovs), length.out = 100)
  p0.posterior1 <- do.call( rbind,
                            lapply(1:length(myResults_M$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_M$sims.list$p0[x,1]) +
                                              myResults_M$sims.list$betaTrap[x] * detCovs)
                                   }))
  p0.posterior2 <- do.call( rbind,
                            lapply(1:length(myResults_M$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_M$sims.list$p0[x,2]) +
                                              myResults_M$sims.list$betaTrap[x] * detCovs)
                                   }))
  p0.posterior3 <- do.call( rbind,
                            lapply(1:length(myResults_M$sims.list$betaTrap),
                                   function(x){
                                     ilogit(logit(myResults_M$sims.list$p0[x,3]) +
                                              myResults_M$sims.list$betaTrap[x] * detCovs)
                                   }))
  
  maxP0 <- round(max(c(p0.posterior1,p0.posterior2,p0.posterior3)),3)    
  p0cols <- hcl.colors(3)
  plot( x = detCovs, y = colMeans(p0.posterior1),
        type = "n", ylim = c(0, maxP0), xlim = range(detCovs),
        ylab = "p0", xlab = "search tracks length (km)", axes = FALSE)
  
  minTracks <- min(trackLength[])
  maxTracks <- max(trackLength[])
  xLabels <- round(seq(minTracks, maxTracks, length.out = 10),2)
  axis(1, at = round(seq(min(detCovs), max(detCovs), length.out = 10),3),
       labels = xLabels, 
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxP0,length.out = 6),
       labels = seq(0,maxP0,length.out = 6),
       tck = 0.01, las = 1)
  
  for(i in seq(1,dim(p0.posterior1)[1],10)){
    points( x = detCovs, y = p0.posterior1[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[1], alpha.f = 0.01))
    points( x = detCovs, y = p0.posterior2[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[2], alpha.f = 0.01))
    points( x = detCovs, y = p0.posterior3[i, ], lwd = 3,
            type = "l", col = adjustcolor(p0cols[3], alpha.f = 0.01))
  }
  
  points( x = detCovs, y = colMeans(p0.posterior1), lwd = 2, 
          type = "l", col = p0cols[1])
  points( x = detCovs, y = colMeans(p0.posterior2), lwd = 2, 
          type = "l", col = p0cols[2])
  points( x = detCovs, y = colMeans(p0.posterior3), lwd = 2, 
          type = "l", col = p0cols[3])
  
  mtext( text = "Males", side = 3, line = -3, font = 2, cex = 1.5)
  
  legend("topright", bty = "n",
         legend = c("BFNP","NR", "SNP"), 
         lwd = 2, col = p0cols)
  # dev.off()
  
  
  
  ## ------   5. PARAMETER TABLES ------
  ## ------     5.1. FEMALES ------
  params.simple <- names(myResults_F$mean)[2:8]
  params.means <- do.call(c, lapply(params.simple, function(x)myResults_F$mean[[x]][1]))
  params.sd <- do.call(c, lapply(params.simple, function(x)myResults_F$sd[[x]][1]))
  params.q2.5 <- do.call(c, lapply(params.simple, function(x)myResults_F$q2.5[[x]][1]))
  params.q97.5 <- do.call(c, lapply(params.simple, function(x)myResults_F$q97.5[[x]][1]))
  # params.Rhat <- do.call(c, lapply(params.simple, function(x)myResults_F$Rhat[[x]][1]))
  # params.n.eff <- do.call(c, lapply(params.simple, function(x)myResults_F$n.eff[[x]][1]))
  
  params.summary <- cbind.data.frame( round(params.means,4)[c(1,2,5,7,3,6,4)],
                                      # round(params.sd,5),
                                      round(params.q2.5,4)[c(1,2,5,7,3,6,4)],
                                      round(params.q97.5,4)[c(1,2,5,7,3,6,4)])#, params.Rhat,params.n.eff)
  
  params.summary[2, ] <- params.summary[2, ] * 4 
  params.summary[3, ] <- params.summary[3, ] * 2 
  params.summary[4, ] <- params.summary[4, ] * 5 
  
  names(params.summary) <- c("mean","2.5%CI", "97.5%CI")
  #"sd", "2.5%CI", "97.5%CI", "Rhat", "n.eff")
  #rownames(params.summary) <- params.simple[c(1,2,5,7,3,6,4)]
  rownames(params.summary) <- c("beta_BFNP", "beta_BFNPp", "beta_NR","beta_SNPp",
                                "beta_elev", "beta_disturb", "beta_interact")
  
  ##-- Export as .pdf
  par(mfrow = c(1,1))
  plot.new()
  grid.table(params.summary)
  mtext(text = "Parameter estimates for females",
        side = 3, outer = T, line = -20, cex = 2, font = 2)
  
  # ##-- Export as .csv
  # write.csv( x = params.summary,
  #            file = file.path(myVars$WD, myVars$outFolder, "TABLES",
  #                             paste( myVars$outFolder, "_params_F.csv", sep = "")))
  # 
  # ##-- Export as .tex
  # print( xtable(params.summary, type = "latex"),
  #        floating = FALSE, # scalebox=.8,
  #        add.to.row = list(list(seq(1, nrow(params.summary), by = 2)),"\\rowcolor[gray]{.95} "),
  #        file = file.path(myVars$WD, myVars$outFolder,"TABLES",
  #                         paste( myVars$outFolder, "_params_F.tex", sep = "")))
  # 
  
  ## ------     5.2. MALES ------
  params.simple <- names(myResults_M$mean)[2:8]
  params.means <- do.call(c, lapply(params.simple, function(x)myResults_M$mean[[x]][1]))
  params.sd <- do.call(c, lapply(params.simple, function(x)myResults_M$sd[[x]][1]))
  params.q2.5 <- do.call(c, lapply(params.simple, function(x)myResults_M$q2.5[[x]][1]))
  params.q97.5 <- do.call(c, lapply(params.simple, function(x)myResults_M$q97.5[[x]][1]))
  # params.Rhat <- do.call(c, lapply(params.simple, function(x)myResults_M$Rhat[[x]][1]))
  # params.n.eff <- do.call(c, lapply(params.simple, function(x)myResults_M$n.eff[[x]][1]))
  
  params.summary <- cbind.data.frame( round(params.means,4)[c(1,2,5,7,3,6,4)],
                                      #round(params.sd,4),
                                      round(params.q2.5,4)[c(1,2,5,7,3,6,4)],
                                      round(params.q97.5,4)[c(1,2,5,7,3,6,4)])#, params.Rhat,params.n.eff)
  
  params.summary[2, ] <- params.summary[2, ] * 4 
  params.summary[3, ] <- params.summary[3, ] * 2 
  params.summary[4, ] <- params.summary[4, ] * 5 
  
  names(params.summary) <- c("mean", "2.5%CI", "97.5%CI")#, "Rhat", "n.eff")
  #rownames(params.summary) <- params.simple[c(1,2,5,7,3,6,4)]
  rownames(params.summary) <- c("beta_BFNP", "beta_BFNPp", "beta_NR","beta_SNPp",
                                "beta_elev", "beta_disturb", "beta_interact")
  ##-- Export as .pdf
  plot.new()
  grid.table(params.summary)
  mtext(text = "Parameter estimates for males",
        side = 3, outer = T, line = -20, cex = 2, font = 2)
  
  ##-- Export as .csv
  # write.csv( x = params.summary,
  #            file = file.path(myVars$WD, myVars$outFolder, "TABLES",
  #                             paste( myVars$outFolder, "_params_F.csv", sep = "")))
  # 
  # ##-- Export as .tex
  # print( xtable(params.summary, type = "latex"),
  #        floating = FALSE,# scalebox=.8,
  #        add.to.row = list(list(seq(1, nrow(params.summary), by = 2)),"\\rowcolor[gray]{.95} "),
  #        file = file.path(myVars$WD, myVars$outFolder,"TABLES",
  #                         paste( myVars$outFolder, "_params_M.tex", sep = "")))
  
  
  graphics.off()
}
##------------------------------------------------------------------------------
