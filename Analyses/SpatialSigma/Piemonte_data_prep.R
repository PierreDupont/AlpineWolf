################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list = ls())


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
library(gridExtra)
library(MetBrewer)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
##-- ANALYSIS DIRECTORY
thisDir = file.path(analysisDir, "SpatialSigma/Piemonte")

##-- HABITAT SPECIFICATIONS
habitat = list( resolution = 20000, 
                buffer = 30000)

##-- DETECTORS SPECIFICATIONS
detectors = list( resolution = 5000,
                  detSubResolution = 1000,
                  samplingMonths = c(10:12,1:4))

##-- NGS DATA SPECIFICATIONS
data = list( sex = c("F","M"),
             status = c("alpha","pup","other")) 

##-- CREATE FOLDER TO STORE PROCESSED DATA
dir.create(file.path(thisDir, "data"), recursive = T)



## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
##---- Polygon of Italy and neighbouring countries
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

##--- Regions 
regions <- read_sf(file.path(dataDir,"GISData/Output_layout/Alpine_Regions.shp"))
regions <- st_transform(x = regions, crs = st_crs(countries))
regions$ID <- as.numeric(as.factor(regions$DEN_UTS))
plot(st_geometry(regions))
regions <- regions[regions$DEN_REG %in% c("Piemonte", "Valle d'Aosta", "Liguria"), ]


##--- Alps
alps <- read_sf(file.path(dataDir,"GISData/Output_layout/Italian_Alps.shp"))
alps <- st_transform(x = alps, crs = st_crs(countries))
plot(st_geometry(alps))

##--- Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_intersection(.,regions) %>%
  st_union() 

##--- SCR grid
SCRGrid <- read_sf(file.path(dataDir,"GISData/SECR_presence_layer_2020_2021/SECR_presence_layer_2021.shp"))
SCRGrid <- st_transform(x = SCRGrid, crs = st_crs(countries))
SCRGrid <- SCRGrid[SCRGrid$`Pres_20-21` == 1, ]


## 
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
plot(st_geometry(temp),col="gray60")
temp <- st_transform(temp, st_crs(studyArea))

plot(st_geometry(countries), col = "gray80", border = F)
plot(st_geometry(temp), col = "gray80",add=T,border = F)
plot(st_geometry(st_intersection(studyArea,countries)),add=T,col="red",border = F)



## ------   2. SEARCH EFFORT DATA ------
##---- Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))

##---- Convert dates
transects$Date <- parse_date_time(transects$date, orders = c('ymd'))
transects$Year <- as.numeric(format(transects$Date,"%Y"))
transects$Month <- as.numeric(format(transects$Date,"%m"))

##---- Plot check
plot(transects, col = "red", add = T)




## ------   3. DNA DATA ------
ngs <- read.csv(file.path(dataDir,"DNA/ngs21032022.csv"))
dim(ngs)

##---- Use lubridate to clean dates
# The orders argument is a character vector containing the possible date-time parsing 
# formats in the order they should be tested. So by giving c('dmy', 'ymd'), lubridate 
# will try to parse all strings as Day, Month, Year format. If it can't do that 
# successfully (for example, the date 2021-03-18 won't work as there is no 2021th day), 
# it will try the next in the list until all strings are parsed, or all orders exhausted.
ngs$Date <- parse_date_time(ngs$Date, orders = c('dmy','mdy','ymd'))
ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
ngs$Month <- as.numeric(format(ngs$Date,"%m"))

##---- Filter out samples from 2019
ngs <- ngs[ngs$Year > 2019, ]
dim(ngs)

##---- Filter out dead recoveries
ngs <- ngs[ngs$Dead.recovery == "", ]
dim(ngs)

##---- Number of detections per individual
numDetsPerId <- table(ngs$Genotype.ID)
hist(numDetsPerId)

##---- Number of individuals detected
length(numDetsPerId)

##---- Mean number of detections per ID
mean(numDetsPerId)

##---- Number of individuals with > 1 detection
sum(numDetsPerId > 1)

##---- mean number of individual detections by sex
table(ngs$Sex, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Sex, useNA = "always"),
      2,
      function(x)sum(x)/sum(x>0))

##---- mean number of individual detections by status
ngs$Status[ngs$Status == "pup 2016"] <- "pup"
ngs$Status[ngs$Status == "pup 2017"] <- "pup"
table(ngs$Status, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Status, useNA = "always"),
      2,
      function(x)sum(x)/sum(x>0))

##---- mean number of individual detections by sex and status
apply(table(ngs$Genotype.ID, ngs$Status, ngs$Sex, useNA = "always"),
      c(2,3),
      function(x)sum(x)/sum(x>0))

##---- Turn ngs into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$CoordX, ngs$CoordY)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

##---- Plot check
plot(st_geometry(studyArea))
plot(st_geometry(countries), col = "gray80",add=T)
plot(st_geometry(transects), col = "red", add = T)
plot(st_geometry(ngs), add = T, pch = 3)



## ------   4. SPATIAL COVARIATES ------
##---- Read snow-fall rasters
tif <- read_stars(file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/ERA-5_snowfall_alps_32N.nc"))
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


##---- Load historical presence data
packPres_files <- list.files(path = file.path(dataDir, "/GISData/Branchi_SECR"),
                             pattern = ".shp")
packPres_raw <- lapply(packPres_files, function(x){
  tmp <- read_sf(file.path(dataDir, "/GISData/Branchi_SECR", x))
  tmp <- st_transform(x = tmp, crs = st_crs(studyArea))
  tmp
})
plot(studyArea)
lapply(packPres_raw, function(x)plot(x, add = T))


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
PA <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Protected_Areas/PA.shp"))
plot(studyArea)
plot(PA, add = T)

## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. DETECTORS ------
## ------     1.1. DETECTORS CHARACTERISTICS ------
##---- Make PAB search grid
searchGrid <- MakeSearchGrid(
  data = as_Spatial(studyArea),
  resolution = detectors$resolution,
  div = (detectors$resolution/detectors$detSubResolution)^2,
  center = T,
  plot = TRUE,
  fasterize = TRUE)

##---- Create detector raster
detectors$raster <- raster(extent(st_bbox(studyArea)))
res(detectors$raster) <- detectors$resolution
detectors$raster <- fasterize(countries, detectors$raster)

##---- Mask and crop raster cells outside the study area to obtain the detector grid
detectors$raster <- mask(detectors$raster, st_as_sf(studyArea))
detectors$raster <- crop(detectors$raster, st_as_sf(studyArea))

##---- Remove glaciers
glaciers <- CLC
glaciers[!glaciers[ ] %in% c(9)] <- 0
glaciers[glaciers[ ] %in% c(9)] <- 1
to.remove <- glaciers %>%
  raster::aggregate( x = .,
                     fact = detectors$resolution/res(glaciers),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>0.5) %>%
  st_as_sf() %>%
  fasterize(.,detectors$raster)

detectors$raster[to.remove[ ] == 1] <- NA
plot(detectors$raster)

##---- Remove big lakes
lakes <- CLC
lakes[!lakes[ ] %in% c(7)] <- 0
lakes[lakes[ ] %in% c(7)] <- 1
to.remove <- lakes %>%
  raster::aggregate( x = .,
                     fact = detectors$resolution/res(lakes),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>0.5) %>%
  st_as_sf() %>%
  fasterize(.,detectors$raster)
detectors$raster[to.remove[ ] == 1] <- NA
plot(detectors$raster)

##---- Remove big cities
to.remove <- POP_raw %>%
  raster::aggregate( x = .,
                     fact = detectors$resolution/res(POP_raw),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>5000) %>%
  st_as_sf() %>%
  fasterize(.,detectors$raster)
detectors$raster[to.remove[ ] == 1] <- NA
plot(detectors$raster)

##---- Transform into spatial grid
detectors$grid <- st_as_sf(rasterToPolygons(detectors$raster))
detectors$grid$id <- 1:nrow(detectors$grid)
st_crs(detectors$grid) <- st_crs(studyArea)

##---- Extract length and number of transects in each grid cell
intersection <- st_intersection(detectors$grid, transects) %>%
  mutate(LEN = st_length(.),
         QI = .$Q.index) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
            transect_N = length(unique(Date)),   ## Get total number of visits in each detector grid cell
            transect_qi = mean(QI))              ## Get mean transects quality index for each detector grid cell

##---- Store in detector grid
detectors$grid <- detectors$grid %>%
  left_join(intersection, by = "id")
detectors$grid$transect_L[is.na(detectors$grid$transect_L)] <- 0
detectors$grid$transect_N[is.na(detectors$grid$transect_N)] <- 1
detectors$grid$transect_qi[is.na(detectors$grid$transect_qi)] <- 0
detectors$grid$transect_L <- scale(detectors$grid$transect_L)
detectors$grid$mean_transect_L <- scale(detectors$grid$transect_L/detectors$grid$transect_N)

# ####---- If you want only grid cells that contain transects instead
# detectors$grid <- detectors$grid %>%
#   left_join(intersection, by = "id") %>%
#  filter(!is.na(transect_L))

detPoly <- rasterToPolygons(detectors$raster)
proj4string(detPoly) <-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"

##---- Remove sub-detectors outside the habitat
detToKeep <- over(x = searchGrid$detector.sp,
                  y = detPoly)
searchGrid$detector.sp <- searchGrid$detector.sp[!is.na(detToKeep[,1]), ]

mainDetToKeep <- searchGrid$main.detector.sp$main.cell.id %in% searchGrid$detector.sp$main.cell.id
searchGrid$main.detector.sp <- searchGrid$main.detector.sp[mainDetToKeep, ]

##---- Re-order detectors to match the detector grid
searchGrid$main.detector.sp <- searchGrid$main.detector.sp[
  order(-searchGrid$main.detector.sp$main.cell.y,
        searchGrid$main.detector.sp$main.cell.x), ]
searchGrid$main.detector.sp$main.cell.new.id <- 1:nrow(detectors$grid)

##---- Re-name main detectors to match the detector grid
searchGrid$detector.sp$main.cell.new.id <- sapply(
  searchGrid$detector.sp$main.cell.id, function(x){
    which(searchGrid$main.detector.sp$main.cell.id == x)
  })

##---- Store detector coordinates
detectors$coords <- st_coordinates(st_centroid(detectors$grid))[ ,c("X","Y")]
dimnames(detectors$coords) <- list(1:nrow(detectors$grid),
                                   c("x","y"))
detectors$sp <- searchGrid$main.detector.sp
detectors$sub.sp <- searchGrid$detector.sp
detectors$grid$subdetector <- as.numeric(table(detectors$sub.sp$main.cell.new.id))

##---- Extract total number of detectors
detectors$n.detectors <- nrow(detectors$grid)




## ------     1.2. DETECTOR COVARIATES ------
## ------       1.2.1. SNOWFALL ERA-5 ------
##---- Calculate mean and cumulative snowfall
SNOW$snow.mean <- rowMeans(st_drop_geometry(SNOW), na.rm = T)
SNOW$snow.sum <- rowSums(st_drop_geometry(SNOW), na.rm = T)

##---- Extract snow-fall in each detector grid cell
intersection <- st_intersection(detectors$grid, SNOW) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(snow_fall = mean(snow.mean))

##---- Store average snow-fall
detectors$grid <- detectors$grid %>%
  left_join(intersection, by = "id")

##---- scale covariate
detectors$grid$snow_fall[is.na(detectors$grid$snow_fall)] <- 0
detectors$grid$snow_fall <- scale(detectors$grid$snow_fall)




## ------       1.2.2. CORINE LAND COVER ------
##---- Identify layers we want to use
layer.index <- c(1,2,3,4,5,9)

##---- Set-up data frame to store CLC values
CLC.df <- as.data.frame(matrix(NA, detectors$n.detectors, length(layer.index)))
names(CLC.df) <- layer.names[layer.index]
CLC.df$id <- 1:detectors$n.detectors
COV <- list()

##---- Loop over the different layers of interest
par(mfrow = c(1,2))
for(l in 1:length(layer.index)){
  temp <- CLC
  temp[!temp[] %in% layer.index[l]] <- 0
  temp[temp[] %in% layer.index[l]] <- 1
  plot(temp)
  plot(st_geometry(detectors$grid), add = T)
  
  COV[[l]] <- raster::aggregate( x = temp,
                                 fact = detectors$resolution/res(temp),
                                 fun = mean)
  plot(COV[[l]])
  plot(st_geometry(detectors$grid), add = T)
  CLC.df[ ,l] <- raster::extract( x = COV[[l]],
                                  y = detectors$sp)
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


## ------     1.3. FORMAT AS DATAFRAME -------
detectors$tdf <- cbind.data.frame("x" = st_coordinates(st_centroid(detectors$grid))[ ,1],
                                  "y" = st_coordinates(st_centroid(detectors$grid))[ ,2],
                       st_drop_geometry(detectors$grid))
head(detectors$tdf)




## ------   2. HABITAT ------
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (detector grid + buffer)
habitat$polygon <- st_buffer(st_union(detectors$grid),
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
plot(habitat$raster)

##---- Remove glaciers
glaciers <- CLC
glaciers[!glaciers[ ] %in% c(9)] <- 0
glaciers[glaciers[ ] %in% c(9)] <- 1
to.remove <- glaciers %>%
  raster::aggregate( x = .,
                     fact = habitat$resolution/res(temp),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>0.5) %>%
  st_as_sf() %>%
  fasterize(.,habitat$raster)

habitat$raster[to.remove[ ] == 1] <- NA
plot(habitat$raster)

##---- Remove big lakes
lakes <- CLC
lakes[!lakes[ ] %in% c(7)] <- 0
lakes[lakes[ ] %in% c(7)] <- 1
to.remove <- lakes %>%
  raster::aggregate( x = .,
                     fact = habitat$resolution/res(lakes),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>0.5) %>%
  st_as_sf() %>%
  fasterize(.,habitat$raster)
habitat$raster[to.remove[ ] == 1] <- NA
plot(habitat$raster)

##---- Remove big cities
to.remove <- POP_raw %>%
  raster::aggregate( x = .,
                     fact = habitat$resolution/res(POP_raw),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>2000) %>%
  st_as_sf() %>%
  fasterize(.,habitat$raster)
habitat$raster[to.remove[ ] == 1] <- NA
plot(habitat$raster)



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
comparison <- st_buffer(st_union(detectors$grid[as.numeric(detectors$grid$transect_L) >0, ]),5000)
habitat$comparison <- mask(habitat$raster, st_as_sf(comparison))
habitat$comparison[habitat$comparison == 0] <- NA

##---- Create "EXTRACTION" habitat raster (for the report)
habitat$extraction <- mask(habitat$raster, st_as_sf(st_union(SCRGrid)))
habitat$extraction[habitat$extraction == 0] <- NA

##---- Create "ITALY" habitat raster (for density estimates)
habitat$Italia <- mask(habitat$raster, st_as_sf(studyArea))
habitat$Italia[habitat$Italia == 0] <- NA

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
plot(transects,add=T, col="red")
plot(st_geometry(detectors$grid), add=T)




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



## ------       2.2.6. HISTORICAL PRESENCE ------
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




## ------       2.2.7. IUCN PRESENCE ------
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




## ------     2.3. FORMAT AS DATAFRAME -------
habitat$ssdf <- cbind.data.frame("x" = st_coordinates(st_centroid(habitat$grid))[ ,1],
                                 "y" = st_coordinates(st_centroid(habitat$grid))[ ,2],
                       st_drop_geometry(habitat$grid))
head(habitat$ssdf)



## ------   3. DETECTION DATA ------
##---- Calculate distance between detections and detectors
ngs <- st_intersection(ngs, studyArea)
closest <- nn2( coordinates(detectors$sub.sp),
                st_coordinates(ngs),
                k = 1,
                searchtype = "radius",
                radius = 10000)


##---- Assign each detection to a detector based on minimum distance
ngs$subdetector <- c(closest$nn.idx)
ngs$detector <- detectors$sub.sp$main.cell.new.id[closest$nn.idx] 

##---- Drop duplicated detections ate the same sub-detectors
ngs <- ngs[!duplicated(ngs[ ,c("subdetector", "Genotype.ID")]), ] 

##---- Count individual detections per detector
edf <- cbind.data.frame( st_coordinates(ngs),
                         ngs$Genotype.ID,
                         ngs$detector,
                         ngs$subdetector)
names(edf) <- c("x","y","id","detector","subdetector")
head(edf)




## ------   4. INDIVIDUAL COVARIATES ------
##---- List detected indivdual names
IDs <- unique(edf$id)

##---- Set-up vectors to store individual covariates
sex <- status <- pack <- rep(NA, length(IDs))
for(i in 1:length(IDs)){
  
  ##-- Sex 
  temp <- unique(ngs$Sex[ngs$Genotype.ID %in% IDs[i]])
  if(length(temp)>1)warning(print(paste("ID:", IDs[i], " sex:", temp)))
  sex[i] <- ifelse(length(temp>1),temp[1],temp)
  
  ##-- Social status
  temp <- unique(ngs$Status[ngs$Genotype.ID %in% IDs[i]])
  if(length(temp) > 1){
    warning(print(paste("ID:", IDs[i], " status:", temp)))
    if(any(temp %in% "pup")){status[i] <- "pup"}
    if(any(temp %in% "other")){status[i] <- "other"}
    if(any(temp %in% "alpha")){status[i] <- "alpha"}
  } else {
    status[i] <- temp
  }
  
  ##-- Pack membership
  temp <- unique(ngs$Pack[ngs$Genotype.ID %in% IDs[i]])
  if(any(temp %in% "dispersal")){
    status[i] <- "other"
  }
}

##---- Convert to numerical values
sex[sex == "F"] <- 0
sex[sex == "M"] <- 1
sex[sex == ""] <- NA
sex <- as.numeric(sex)

status[status == "alpha"] <- 1
status[status == "pup"] <- 2
status[status == "na"] <- NA
status[status == "other"] <- 3

status <- as.numeric(status)

id.covs <- cbind.data.frame( "id" = IDs,
                             "sex" = sex,
                             "status" = status)
head(id.covs)





## ------   5. SAVE PREPARED DATASET -----
save( habitat, detectors, ngs, edf, id.covs,
      file = file.path(thisDir, "data", "Piemonte_data.RData"))



