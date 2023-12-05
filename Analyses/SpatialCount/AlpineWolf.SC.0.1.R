## -------------------------------------------------------------------------- ##
## ----------------------- ALPINE WOLF SCR ---------------------------------- ##
## ----------------- s[CLC + HumPop + Zone] --------------------------------- ##
## ---------------- z[psi,rho,theta] ---------------------------------------- ##
## ---------------- y[p0(transects + zone + ),sigma] ------------------------ ##
## -------------------------------------------------------------------------- ##
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
library(stars)
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
library(MCMCvis)



## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SC.0.1"
thisDir <- file.path(analysisDir, modelName)


if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}


## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 30000)

# ## DETECTORS SPECIFICATIONS
# detectors = list( resolution = 5000,
#                   detSubResolution = 1000,
#                   samplingMonths = c(10:12,1:4))

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
grid.r <- st_as_stars(st_bbox(studyArea), dx = habitat$resolution, dy = habitat$resolution)
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
  st_geometry() # since you want the centroids in a second geometry col
# plot(studyAreaGrid["SCR"], add = T)

count <- read_sf(file.path(dataDir,
                           "GISData/Countries/Countries_WGS84.shp"))
temp <- count[count$CNTRY_NAME %in% c("Albania",
                                      "Austria",
                                      "Bosnia and Herzegovina",
                                      "Bulgaria",
                                      #"Czech Republic",
                                      "France",
                                      "Germany",
                                      "Greece",
                                      "Croatia",
                                      "Hungary",
                                      "Macedonia",
                                      "Malta",
                                      "Montenegro",
                                      "Serbia",
                                      #"Romania",
                                      "Slovenia",
                                      "San Marino",
                                      "Switzerland"), ]
# plot(st_geometry(temp),col="gray60")
temp <- st_transform(temp, st_crs(studyArea))


## ------   2. CAMERA TRAPS DATA ------
##---- Load GPS of Camera Traps
ct <- read_sf(file.path(dataDir,"GISData/CameraTraps/Ctraps/ct_alpi_231201.shp"))

##---- Convert dates
ct$date_st <- parse_date_time(ct$date_st, orders = c('dmy'))
ct$date_en <- parse_date_time(ct$date_en, orders = c('dmy'))


ct$year_st <- as.numeric(format(ct$date_st,"%y"))
ct$month_st <- as.numeric(format(ct$date_st,"%m"))
ct$year_en <- as.numeric(format(ct$date_en,"%y"))
ct$month_en <- as.numeric(format(ct$date_en,"%m"))
ct$tot_attivi <- as.numeric(ct$tot_attivi) 

##---- Plot check
# plot(studyArea, col="steelblue")
# plot(ct$geometry, col = "red", pch=16, add=T)

##---- Number of pics per CT
hist(ct$tot_attivi)

ct_cl <- ct[,c(1,6)]
ct_cl$id <- paste0("FT", 1:nrow(ct_cl))

ct_cl$coords <- st_coordinates(ct_cl)
colnames(ct_cl$coords) = c("x", "y")


dimnames(ct_cl$coords) <- list(1:nrow(ct_cl$coords),
                                   c("x","y"))

## ------   3. PICTURES/SIGHTENING DATA ------
##---- Images from Camera Traps data
pics_raw <- read_sf(file.path(dataDir,"GISData/CameraTraps/photos/ft_alpi_photos_231201.shp"))
dim(pics_raw)
# plot(studyArea, col="steelblue")
# plot(pics_raw$geometry,  col = "blue", pch=16, add=T)


##---- Use lubridate to clean dates
# The orders argument is a character vector containing the possible date-time parsing 
# formats in the order they should be tested. So by giving c('dmy', 'ymd'), lubridate 
# will try to parse all strings as Day, Month, Year format. If it can't do that 
# successfully (for example, the date 2021-03-18 won't work as there is no 2021th day), 
# it will try the next in the list until all strings are parsed, or all orders exhausted.
pics_raw$date <- parse_date_time(pics_raw$date, orders = c('dmy','mdy','ymd'))
pics_raw$Year <- as.numeric(format(pics_raw$date,"%Y"))
pics_raw$Month <- as.numeric(format(pics_raw$date,"%m"))
pics_raw$n_lupi <- as.numeric(pics_raw$n_wolves)

##---- Filter out samples outside the monitoring period
pics_1 <- filter(pics_raw,Month %in% c(1,2,3,4,10,11,12)) 
pics_na <- pics_raw[is.na(pics_raw$Month),]
pics_t <- rbind(pics_1,pics_na)
dim(pics_t)
pics_t <- filter(pics_t,C1 %in% 'C1')  
dim(pics_t)

pics_t <- st_join(ct_cl, pics_t)


##---- Aggregate pictures per DAY ----
pics_d <- pics_t %>% 
  group_by(day = lubridate::day(date), id) %>% 
  summarise(n_wolves = max(n_wolves))
# hist(pics_d$n_wolves)
dim(pics_d)


##---- Aggregate pictures per WEEK ----
pics_w <- pics_t %>% 
  group_by(week = lubridate::week(date), id) %>% 
  summarise(n_wolves = max(n_wolves))
# hist(pics_w$n_wolves)
dim(pics_w)


##---- Aggregate pictures per MONTH ----
pics_m <- pics_t %>% 
  group_by(Month, id) %>% 
  summarise(n_wolves = max(n_wolves))
# hist(pics_m$n_wolves)
dim(pics_m)



##---- Number of pics per CT -----
# pick whatever dataset you want to use between the aggregations above (months, week, day, raw)
pics <- pics_t

# same info -----
numDetsPerCT <- table(pics$id)
# which(is.na(as.numeric(pics_t$id_ct)))
hist(numDetsPerCT)

##---- Number of at least one detection
length(numDetsPerCT)

##---- Mean number of detections per CT
mean(numDetsPerCT)

##---- Number of CT with > 1 detection
sum(numDetsPerCT > 1)
max(numDetsPerCT)


##---- Maximum number of individuals per CT
maxDetsPerCT <- pics %>% group_by(id) %>% 
  slice_max(n_wolves, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(id, n_wolves) %>% 
  group_by(id) %>% 
  summarise(n_wolves = sum(n_wolves))
# Plot
hist(maxDetsPerCT$n_wolves)


#plot(maxDetsPerCT)

##---- mean number of individual detections by CT
MeanDet <- apply(table(pics$n_wolves, pics$id, useNA = "always"), 2, function(x)sum(x)/sum(x>0))
hist(MeanDet)

##---- DATA WRANGLING -----
# Filter out pics with no camera trap ID
pics_l <- subset(pics, trimws(id) !="")

# Keep only ct_it, and number of wolves detected
pics_l <- pics_l[, c(2,4,19)]
# Give index based on detections per CT
pics_l2 <- pics_l %>%
  group_by(id) %>%
  mutate(index = row_number()) %>%
  ungroup()

det_w <- pics_l2 %>%
  pivot_wider(names_from = "index", values_from = "n_wolves")


det_w <- det_w %>%
  mutate(tot_wolves = rowSums(pick(4:68), na.rm = TRUE),
         tot_detections = rowSums(!is.na(pick(4:68))))

det_w <- det_w[,c(2,1,69:70,4:68,3)]

##---- Plot check
# plot(st_geometry(studyArea), col="steelblue")
# # plot(st_geometry(countries), col = "gray80",add=T)
# plot(st_geometry(ct), col = "red", add = T)
# plot(st_geometry(pics), add = T, pch = 3)

## ------   4. SPATIAL COVARIATES ------
##---- Read snow-fall rasters
tif <- read_stars(file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/ERA-5_snowfall_alps_32N.tif"))
SNOW <- st_as_sf(tif)


##---- Load Corine Land Cover data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))
## 1 = developed
## 2 = agriculture
## 3 = forest
## 4 = herbaceous
## 5 = bare rock
## 6 = sea features
## 7 = inland water
## 9 = perpetual snow
layer.names <- c("developed","agriculture","forest","herbaceous","bare rock",
                 "sea features","inland water","NA","perpetual snow")
plot(CLC)


##---- Load elevation data
DEM_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))
plot(DEM_raw)


##---- Load terrain ruggedness data
TRI_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))
plot(TRI_raw)


##---- Load human population density data
POP_raw  <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
plot(POP_raw)


##---- Load roads data
roads <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Road Density/road_density_tot_clipped.shp"))
table(roads$highway) # see https://wiki.openstreetmap.org/wiki/Key:highway for road classes description
mainRoads <- subset(roads, highway %in% c("motorway", "trunk", "primary"))
rm(roads)


##---- Create a line to delineate the "western Alps"
cutline <- st_linestring(rbind(c(700000,4800000), c(500000,5230000)))
cutline <- st_sfc(cutline)
st_crs(cutline) <- st_crs(studyArea)
cutline <- st_buffer(cutline, dist = 0.0001)
plot(studyArea)
plot(cutline, add = T, border = "blue", lwd = 2)

##---- Create east and western Alps polygons
studyArea_west <- st_sfc(st_cast(st_difference(st_buffer(studyArea,50000), cutline),"POLYGON")[[1]])
st_crs(studyArea_west) <- st_crs(studyArea)
studyArea_east <- st_difference(st_buffer(studyArea,50000),studyArea_west)


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


##---- Load protected areas
# PA <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Protected_Areas/PA.shp"))
# plot(studyArea)
# plot(PA, add = T)



## ------   5. PRE-PROCESSED STUFF ------
load(file.path(thisDir,"Habitat.RData"))
load(file.path(thisDir,"Detectors.RData"))
 

## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. DETECTORS ------
## ------     1.1. DETECTORS CHARACTERISTICS ------
##---- Make PAB search grid
# Buffer around camera traps
# ct_cl$buffer <- st_buffer(ct_cl, dist = 50)
# Define detectors
detectors <- ct_cl[,c(4,2,5)]


# plot(st_geometry(studyArea), col="steelblue")
# # plot(st_geometry(countries), col = "gray80",add=T)
# plot(detectors$geometry, col="red", add=T)
# plot(st_geometry(pics), add = T, pch = 3)


detectors$n.detectors <- 1:nrow(detectors)
# st_crs(detectors$buffer) <- st_crs(studyArea)



## ------     1.2. DETECTOR COVARIATES ------
## ------       1.2.1. SNOWFALL ERA-5 ------
##---- Calculate mean and cumulative snowfall
SNOW$snow.mean <- rowMeans(st_drop_geometry(SNOW), na.rm = T)
SNOW$snow.sum <- rowSums(st_drop_geometry(SNOW), na.rm = T)

##---- Extract snow-fall in each detector grid cell
intersection <- st_intersection(detectors$buffer, SNOW) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(snow_fall = mean(snow.mean))

##---- Store average snow-fall
detectors$grid <- detectors$buffer %>%
  left_join(intersection, by = "id")

##---- scale covariate
detectors$buffer$snow_fall[is.na(detectors$buffer$snow_fall)] <- 0
detectors$buffer$snow_fall <- scale(detectors$buffer$snow_fall)




## ------       1.2.2. CORINE LAND COVER ------
##---- Identify layers we want to use
layer.index <- c(1,2,3,4,5,9)

##---- Set-up data frame to store CLC values
CLC.df <- as.data.frame(matrix(NA, detectors$n.detectors, length(layer.index)))
names(CLC.df) <- layer.names[layer.index]
CLC.df$id <- 1:detectors$n.detectors
COV <- list()

##---- Loop over the different layers of interest
for(l in 1:length(layer.index)){
  temp <- st_as_stars(CLC)
  temp[!temp[] %in% layer.index[l]] <- 0
  temp[temp[] %in% layer.index[l]] <- 1
  # plot(temp)
  # plot(st_geometry(detectors$geometry), add = T)

  # COV[[l]] <- raster::aggregate(x = temp, detectors$grid)
  # plot(COV[[l]])
  # plot(st_geometry(detectors$grid), add = T)
  CLC.df[ ,l] <- st_extract( x = temp,
                             at = detectors$geometry)
  print(layer.index[l])
}#c

##---- Store in detector grid
detectors$grid <- left_join(detectors$grid, CLC.df, by = "id")


## ------       1.2.3. ELEVATION -------
##---- Aggregate to the detector resolution
DEM <- raster::aggregate( x = DEM_raw,
                          fact = detectors$resolution/res(DEM_raw),
                          fun = mean)
plot(DEM)
plot(st_geometry(detectors$grid), add = T)

##---- Extract scaled elevation values
detectors$grid$elev <- DEM %>%
  raster::extract(y = detectors$sp) %>%
  scale()




## ------       1.2.4. TERRAIN RUGGEDNESS INDEX -------
##---- Aggregate to the detector resolution
TRI <- raster::aggregate( x = TRI_raw ,
                          fact = detectors$resolution/res(TRI_raw ),
                          fun = mean)

plot(TRI)
plot(st_geometry(detectors$grid), add = T)

##---- Extract scaled terrain ruggedness
detectors$grid$tri <- TRI %>%
  raster::extract(y = detectors$sp) %>%
  scale()




## ------       1.2.5. HUMAN POPULATION DENSITY -------
##---- Aggregate to the detector resolution
POP <- raster::aggregate( x = POP_raw ,
                          fact = detectors$resolution/res(POP_raw ),
                          fun = mean)

plot(POP)
plot(st_geometry(detectors$grid), add = T)

##----Extract scaled human pop density
detectors$grid$pop <- POP %>%
  raster::extract(y = detectors$sp) %>%
  scale()

##---- Extract log(human pop density + 1)
log_POP <- POP
log_POP[ ] <- log(POP[]+1)
detectors$grid$log_pop <- log_POP %>%
  raster::extract(y = detectors$sp) %>%
  scale()




## ------       1.2.6. ROAD DENSITY -------
##---- Extract road length in each detector grid cell
intersection <- st_intersection(detectors$grid, mainRoads) %>%
  mutate(length = st_length(.))  %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(mainRoads_L = sum(length))

##---- Store scaled road density
detectors$grid <- detectors$grid %>%
  left_join(intersection, by = "id")
detectors$grid$mainRoads_L[is.na(detectors$grid$mainRoads_L)] <- 0
detectors$grid$mainRoads_L <- scale(detectors$grid$mainRoads_L)


##---- Display all detector covariates
plot(detectors$grid[ ,3:length(detectors$grid)], max.plot = 15)




## ------       1.2.7. EAST/WEST ------
##---- Calculate distance of all detectors to each zone
distWest <- st_distance(detectors$grid,studyArea_west)
distEast <- st_distance(detectors$grid,studyArea_east)

##---- Assign detectors to the closest zone
detectors$grid$zone <- apply(cbind(distWest,distEast),1,
                             function(x)which.min(x)-1)




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

##---- Remove glaciers
# glaciers <- CLC
# glaciers[!glaciers[ ] %in% c(9)] <- 0
# glaciers[glaciers[ ] %in% c(9)] <- 1
# to.remove <- glaciers %>%
#   raster::aggregate( x = .,
#                      fact = habitat$resolution/res(temp),
#                      fun = mean) %>%
#   rasterToPolygons(.,fun = function(x)x>0.5) %>%
#   st_as_sf() %>%
#   fasterize(.,habitat$raster)
# 
# habitat$raster[to.remove[ ] == 1] <- NA
# plot(habitat$raster)
# 
# ##---- Remove big lakes
# lakes <- CLC
# lakes[!lakes[ ] %in% c(7)] <- 0
# lakes[lakes[ ] %in% c(7)] <- 1
# to.remove <- lakes %>%
#   raster::aggregate( x = .,
#                      fact = habitat$resolution/res(lakes),
#                      fun = mean) %>%
#   rasterToPolygons(.,fun = function(x)x>0.5) %>%
#   st_as_sf() %>%
#   fasterize(.,habitat$raster)
# habitat$raster[to.remove[ ] == 1] <- NA
# plot(habitat$raster)
# 
# ##---- Remove big cities
# to.remove <- POP_raw %>%
#   raster::aggregate( x = .,
#                      fact = habitat$resolution/res(POP_raw),
#                      fun = mean) %>%
#   rasterToPolygons(.,fun = function(x)x>5000) %>%
#   st_as_sf() %>%
#   fasterize(.,habitat$raster)
# habitat$raster[to.remove[ ] == 1] <- NA
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

##---- Create "COMPARISON" habitat raster (for comparison between models)
# comparison <- st_buffer(st_union(detectors$grid[as.numeric(detectors$grid$transect_L) >0, ]),5000)
# habitat$comparison <- mask(habitat$raster, st_as_sf(comparison))
# habitat$comparison[habitat$comparison == 0] <- NA
# 
# ##---- Create "EXTRACTION" habitat raster (for the report)
# habitat$extraction <- mask(habitat$raster, st_as_sf(st_union(SCRGrid)))
# habitat$extraction[habitat$extraction == 0] <- NA
# 
# ##---- Create "ITALY" habitat raster (for density estimates)
# habitat$Italia <- mask(habitat$raster, st_as_sf(studyArea))
# habitat$Italia[habitat$Italia == 0] <- NA

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
# plot(habitat$raster)
# plot(st_geometry(countries), add = T)
# plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))
# plot(ct_cl,add=T, col="red")



## ------     2.2. HABITAT COVARIATES ------
## ------       2.2.1. CORINE LAND COVER ------
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

habitat$grid$alpine <- habitat$grid$herbaceous + habitat$grid$`bare rock`
habitat$grid$human <- habitat$grid$developed + habitat$grid$agriculture


## ------       2.2.2. ELEVATION -------
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



## ------       2.2.6. EAST/WEST ------
##---- Calculate distance of all detectors to each zone
distWest <- st_distance(habitat$grid, studyArea_west)
distEast <- st_distance(habitat$grid, studyArea_east)

##---- Assign detectors to the closest zone
habitat$grid$zone <- apply(cbind(distWest,distEast),1,
                           function(x)which.min(x)-1)


## ------       2.2.7. HISTORICAL PRESENCE ------
habitat$grid$presence <- 0

for(y in 1:length(packPres_raw)){
  ## Extract historical pack presence in each habitat grid cell
  intersection <- st_intersection(habitat$grid, packPres_raw[[y]]) %>%
    mutate(pres = 1)  %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(sumpres = sum(pres))

  tmp <- habitat$grid %>%
    left_join(intersection, by = "id")
  tmp$sumpres[is.na(tmp$sumpres)] <- 0
  habitat$grid$presence <- habitat$grid$presence + tmp$sumpres
  plot(habitat$grid[,"presence"])
}

habitat$grid$presence <- scale(habitat$grid$presence)




## ------       2.2.8. IUCN PRESENCE ------
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

##---- Identify local detectors
localObjects <- getLocalObjects(
  habitatMask = habitat$binary,
  coords = scaledCoords$coordsDataScaled[ ,1:2],
  dmax = 25)



## ------   4. DETECTION DATA ------
## ------     4.1. DETECTION MATRIX : y ------
##---- Calculate distance between detections and sub-detectors
# closest <- nn2( coordinates(detectors$coords),
#                 st_coordinates(pics),
#                 k = 1,
#                 searchtype = "radius",
#                 radius = 10000)
# 
# ##---- Assign each detection to a main detector based on minimum distance
# pics$sub.detector <- c(closest$nn.idx)
# pics$detector <- detectors$sub.sp$main.cell.new.id[closest$nn.idx] 
# 
# ##---- Drop duplicated detections ate the same sub-detectors
# # pics <- pics[!duplicated(pics[ ,c("sub.detector", "id_ct.x")]), ] 
# #ngs <- droplevels(ngs)
# 
# ##---- Count individual detections per detector
# detMat <- as.matrix(table(pics$id_ct.x, pics$sub.detector))
# 
# # ##---- Convert to binary detections (might keep later for binomial model)
# # detMat[detMat > 0] <- 1
# 
# ##---- Retrieve the number of detectors at which each individual is detected
# detNums <- apply(detMat, 1, function(x) sum(x>0))
# 
# ##--- Set-up matrices to store individual detection frequencies and indices
# n.detected <- dim(detMat)[1]
# detIndices <- matrix(-1, n.detected, max(detNums)*2)
# ySparse <- matrix(-1, n.detected, max(detNums)*2)
# 
# ##---- Fill in the matrix
# for(i in 1:n.detected){
#   ##-- Get where detections occur (detectors index)
#   detIndices[i,1:detNums[i]] <- as.numeric(names(which(detMat[i, ] > 0)))
#   ##-- Get detection frequencies (always 1 in case of Bernoulli) 
#   ySparse[i,1:detNums[i]] <- detMat[i,which(detMat[i, ] > 0)]
# }
# 
# ##---- Combine individual detection number, frequencies and indices
# yCombined <- cbind(detNums, ySparse, detIndices)


## -----------------------------------------------------------------------------
## ------ IV. NIMBLE ------- 
## ------   1. MODEL ------
modelCode1 <- nimbleCode( {

    sigma ~ dunif(0,10)#dgamma(24,8) # weakly informative prior - 38-619 km2
    lam0 ~ dunif(0,10)
    psi ~ dunif(0,1)#dbeta(1,1)
    
    ## Intensity of the AC distribution point process
    # habIntensity[1:numHabWindows] <- exp(habCoeffSlope * habCovs[1:numHabWindows])
    sumHabIntensity <- sum(habIntensity[1:n.habWindows])
    logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
    logSumHabIntensity <- log(sumHabIntensity)
    
    # loop over groups
    for(g in 1:G) {
      
      z[g] ~ dbern(psi)
      
      
      s[g,1:2] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
        upperCoords = upperHabCoords[1:n.habWindows,1:2],
        logIntensities = logHabIntensity[1:n.habWindows],
        logSumIntensity = logSumHabIntensity,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
      
      for(j in 1:J) {
        
        distsq[g,j] <- (s[g,1] - trapCoords[j,1])^2 + (s[g,2] - trapCoords[j,2])^2
        lam[g,j] <- lam0 * exp(-distsq[g,j] / (2*sigma^2)) * z[g] 
      }
    } # End of 1:M
    
    lambda.oper ~ dgamma(1,1)
    
    for(j in 1:J) {
      
      oper[j] ~ T(dpois(lambda.oper),1,)
        bigLambda[j] <- sum(lam[1:G,j]) 
        n[j] ~ dpois(bigLambda[j]*oper[j])
      }
    
    N <- sum(z[1:G])
    D <- N/area
})



## ------   2. BUNDLE DATA ------
##---- Set model constants, data & parameter simulated values (==inits)
G <- 2000
area <- st_area(grid) %>%
  drop_units() %>%
  sum()


nimInits <- list(sigma = rnorm(1,5),
                 lam0 = runif(1),
                 lambda.oper = 1,
                 z = rbinom(G,1,0.6),
                 psi = 0.6)

nimData <- list(n = det_w$tot_wolves,
                habIntensity = rep(1,habitat$n.HabWindows),
                area = area,
                oper = det_w$tot_attivi,
                habitatGrid = localObjects$habitatGrid,
                lowerHabCoords = habitat$loScaledCoords, 
                upperHabCoords = habitat$upScaledCoords)

nimConstants <- list(G = G,
                     J = nrow(det_w),
                     trapCoords = detectors$scaledCoords,
                     # K=ncol(det_w),
                     n.habWindows = habitat$n.HabWindows,
                     y.max = dim(habitat$matrix)[1],
                     x.max = dim(habitat$matrix)[2]) 


nimParams <- c("N", "D", "lam0", "sigma","psi","z","s")


##---- Create the nimble model object
nimModel <- nimbleModel( code = modelCode1,
                         constants = nimConstants,
                         inits = nimInits,
                         data = nimData,
                         check = FALSE,
                         calculate = FALSE) 
nimModel$simulate("s")
nimModel$simulate("oper")

nimModel$calculate()
Cmodel <- compileNimble(nimModel)
modelConf <- configureMCMC(nimModel,
                           thin = 1)

modelConf$addMonitors(nimParams)
modelMCMC <- buildMCMC(modelConf)

CmodelMCMC <- compileNimble(modelMCMC, project = nimModel)
out2 <- runMCMC(CmodelMCMC, 
                niter = 20000,
                nchains = 3)

#   return(as.mcmc(out1))
# })

str(out1)
### check on the first batch of results, 200% likely need to keep running ####
# out.mcmc <- as.mcmc(out1)

MCMCplot(object = out1, 
         params = 'N')


MCMCtrace(object = out1,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "N")

# ##---- Compile the nimble model object to C++
# CsimModel <- compileNimble(nimModel)
# CsimModel$calculate()
# 
# ##---- Configure and compile the MCMC object 
# conf <- configureMCMC( model = nimModel,
#                        monitors = nimParams,
#                        thin = 1)
# 
# ##---- Configure reversible jump
# configureRJ( conf =  conf,
#             targetNodes = 'betaHab.raw',
#             indicatorNodes = 'zRJ',
#             control = list(mean = 0, scale = .2))
# 
# Rmcmc <- buildMCMC(conf)
# compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
# Cmodel <- compiledList$model
# Cmcmc <- compiledList$mcmc
# 
# ##---- Run nimble MCMC in multiple bites
# for(c in 1:1){
#   print(system.time(
#     runMCMCbites( mcmc = Cmcmc,
#                   bite.size = 50,
#                   bite.number = 2,
#                   path = file.path(thisDir, paste0("output/chain",c)))
#   ))
# }
# 




## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 10)

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput$samples)
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



# ## ------     2.2. DENSITY EFFECT PLOT ------
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
# ## EE
# pdf(file = file.path(thisDir, paste0(modelName,"_ridge.pdf")),
#     width = 30, height = 22)
# ridgeMap( ital.R,
#           line.col = "white",
#           fill.col = "lightblue4",
#           grid.col = "gray80",
#           scale = 8,
#           lwd = 2,
#           plot.margin = unit(c(10,10,20,0),"pt"),
#           caption = expression(paste(bold("Alpine Wolf ("),
#                                      italic("Canis lupus"),
#                                      bold(")"))),
#           plot.caption = element_text( face = "bold.italic",
#                                        size = 50,
#                                        colour = "white",
#                                        hjust = 0.02))
# graphics.off()

# g = st_graticule(countries)
# plot(st_geometry(g), axes = TRUE)

##------------------------------------------------------------------------------


##---- PROCESS THE OUTPUT
n.chains <- length(nimOutput$samples)
n.iterations <- dim(nimOutput$samples[[1]])[1]
covNames <- colnames(nimData$hab.covs)
betas.wide <- data.table(res$sims.list$betaHab)
dimnames(betas.wide) <- list(NULL, covNames)

zRJ.wide <- data.table(res$sims.list$zRJ)
dimnames(zRJ.wide) <- list(NULL, covNames)

mods <- apply(zRJ.wide, 1, function(x) {
  paste(covNames[x == 1], collapse = "+")
})
betas.wide$model <- zRJ.wide$model <-
  gsub("\\(Intercept\\)\\+", "", mods)


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

betas.df <- betas.df[order(betas.df$variable, betas.df$model, betas.df$chain),]
zRJ.df <- zRJ.df[order(zRJ.df$variable, zRJ.df$model, zRJ.df$chain),]

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
  color = factor(chain)
)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  xlab("Iteration") +  theme(legend.position = "none")




##---- COEFFICIENT TRACE PLOTS (MODEL-SPECIFIC)
ggplot(data = betas.df, aes(
  x = iteration.model,
  y = value,
  color = factor(chain)
)) +
  geom_line() +
  facet_grid(variable ~ model, margins = FALSE, scales = "free") +
  xlab("Iteration") +  theme(legend.position = "none")




##---- PLOT COEFFICIENT ESTIMATES (OVERALL)
ggplot(betas.df, aes(value, variable, alpha = p.inclusion)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "turquoise",
    color = grey(1)
  ) +
  geom_vline(xintercept = 0)




##---- PLOT COEFFICIENT ESTIMATES (MODEL-SPECIFIC)
ggplot(betas.df, aes(value, variable, alpha = weight)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "magenta",
    color = grey(1)
  ) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ model)
