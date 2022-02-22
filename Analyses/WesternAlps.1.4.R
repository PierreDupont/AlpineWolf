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
modelName = "WesternAlps.1.4"
thisDir <- file.path(analysisDir, modelName)

## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 5000)

## DETECTORS SPECIFICATIONS
detectors = list( detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

## NGS DATA SPECIFICATIONS
data = list( sex = c("female","male"),
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
studyArea <- st_union(studyAreaGrid)

## Plot check
plot(studyArea, col = "gray80")
plot(countries["CNTR_CODE"], col = "gray80", add = T)
plot(studyAreaGrid["SCR"], add = T)




## ------   2. SEARCH EFFORT DATA ------
## Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))

## Convert dates
transects$Date <- as.POSIXct(strptime(transects$date, "%Y-%m-%d"))
transects$Year <- as.numeric(format(transects$Date,"%Y"))
transects$Month <- as.numeric(format(transects$Date,"%m"))

## Plot check
plot(transects, col = "red", add = T)




## ------   3. DNA DATA ------
# allSamples <- read_sf(file.path(dataDir,"GISData/scats/merged_scats.shp"))
# plot(allSamples, add = T)
ngs <- read.csv(file.path(dataDir,"DNA/wa_genetic_samples_Italian_Alps_16_02_2022.csv"))
dim(ngs)

## Use lubridate to clean dates
# The orders argument is a character vector containing the possible date-time parsing 
# formats in the order they should be tested. So by giving c('dmy', 'ymd'), lubridate 
# will try to parse all strings as Day, Month, Year format. If it can't do that 
# successfully (for example, the date 2021-03-18 won't work as there is no 2021th day), 
# it will try the next in the list until all strings are parsed, or all orders exhausted.
ngs$Date <- parse_date_time(ngs$Date, orders = c('dmy', 'ymd'))
ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
ngs$Month <- as.numeric(format(ngs$Date,"%m"))

## Filter out samples from 2019
ngs <- ngs[ngs$Year > 2019, ]
dim(ngs)

## Filter out dead recoveries
ngs <- ngs[ngs$Dead.recovery != "Yes", ]
dim(ngs)

## Explore data
numDetsPerId <- table(ngs$Genotype.ID)

## Number of individuals detected
length(numDetsPerId)

## Mean number of detections per ID
mean(numDetsPerId)

## Number of individuals with > 1 detection
sum(numDetsPerId > 1)

## Number of detections per pack
table(ngs$Pack, useNA = "always")

## mean number of detections per sex
table(ngs$Sex, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Sex, useNA = "always"), 2, function(x)sum(x)/sum(x>0))

## mean number of detections per status
ngs$Status[ngs$Status == "pup 2016"] <- "pup"
ngs$Status[ngs$Status == "pup 2017"] <- "pup"
table(ngs$Status, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Status,useNA = "always"), 2, function(x)sum(x)/sum(x>0))


## Turn into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$Coordinate.WGS84.UTM.East,
                                     ngs$Coordinate.WGS84.UTM.North)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

## Plot check
plot(st_geometry(ngs))
plot(studyArea, add = T)
plot(st_geometry(transects), col = "red", add = T)
plot(st_geometry(ngs), add = T, pch = 3)

# ## Calculate distance to the closest transect
# ## For later: can be quite long with many points and transects
# dist <- st_distance(ngs,transects)
# ngs$dist_to_transect <- apply(dist,1,min)
# range(ngs$dist_to_transect)




## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. DETECTORS ------
## For now: detector grid based on transect locations with the chosen resolution
## Snow fall, transect length an number of visits, is extracted for each detector cell
## each month. Only detector grid cells with counts > 0 are kept (i.e. where people logged GPS tracks)

## Alternative for faster model:
## Collapse all months and use average snowfall, total search length and total number of visits


## ------     1.1. DETECTORS CHARACTERISTICS ------
## Create a line to delineate the "western Alps"
cutline <- st_linestring(rbind(c(630000,4800000),
                               c(430000,5230000)))
cutline <- st_sfc(cutline)
st_crs(cutline) <- st_crs(studyArea)
cutline <- st_buffer(cutline, dist = 0.0001)
plot(studyArea)
plot(cutline,add=T, border = "red", lwd = 2)

##---- Create study area polygon for the western Alps (study area cut to the line) 
studyArea_west <- st_sfc(st_cast(st_difference(st_buffer(studyArea,10000), cutline),"POLYGON")[[1]])
st_crs(studyArea_west) <- st_crs(studyArea)
plot(studyArea_west, add=T, col = "gray60")

##---- Create a detector grid from transects with the desired resolution
transects_west <- st_intersection( transects, studyArea_west)

grid <- st_as_sf(st_make_grid( transects_west,
                               detectors$detResolution))
grid$id <- 1:nrow(grid)


##---- Extract length and number of transects in each grid cell
intersection <- st_intersection(grid, transects_west) %>%
  mutate(LEN = st_length(.))  %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
            transect_N = length(unique(Date)))   ## Get total number of visits in each detector grid cell

##---- Store in detector grid
grid <- grid %>% 
  left_join(intersection, by = "id") %>% 
  filter(!is.na(transect_L))
detectors$grid <- grid

detectors$grid$mean_transect_L <- scale(detectors$grid$transect_L/detectors$grid$transect_N)

##---- Extract detector coordinates
detectors$coords <- st_coordinates(st_centroid(detectors$grid))[ ,c("X","Y")]
dimnames(detectors$coords) <- list(1:nrow(detectors$grid),
                                   c("x","y"))

detectors$trials <- detectors$grid$transect_N

##---- Extract total number of detectors
detectors$n.detectors <- nrow(detectors$grid)
plot(detectors$grid)




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

plot(detectors$grid)




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
comparison <- st_buffer(st_union(detectors$grid),5000) 
habitat$comparison <- mask(habitat$raster, st_as_sf(comparison))
habitat$comparison[habitat$comparison == 0] <- NA


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
plot(st_geometry(countries), col = NA, add = T)
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
  sex[i] <- unique(ngs$Sex[ngs$Genotype.ID %in% IDs[i]])
  
  ## Social status
  status[i] <- unique(ngs$Status[ngs$Genotype.ID %in% IDs[i]])
  # 
  # ## Pack membership
  # pack[i] <- unique(ngs$Pack[ngs$Genotype.ID %in% IDs[i]])
  # 
  # ## ...
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


## ------     2.5. ROAD DENSITY 
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
plot( transects_west, col = "firebrick3", add = T)
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
plot(transects_west, col = "gray80", add = T)
plot(centroids,
     add = T, pch = 21, cex = 1,
     bg = adjustcolor("red",0.2),
     col = "red")

mtext( text = "Individual centroids", side = 3, font = 2)


##-- NGS samples per sex
par(mfrow = c(1,2))
sexCol <- c(hcl.colors(2),"gray60")
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
plot(transects_west, col = "black", add = T)
plot(st_geometry(ngs),
     add = T, pch = 21, cex = 1,
     bg = ifelse(ngs$Sex == "F", sexCol[1], sexCol[2]),
     col = "black")
legend(  x = 502000, y = 5110000, bty = "n",
         legend = c("Females", "males", "NA"),
         title = "NGS samples",
         pch = 21,
         pt.bg = sexCol)

barplot( cbind(table(sex,useNA = "always")),
         col = sexCol,
         ylab = "Number of individuals detected",
         beside = F, width = 1, xlim = c(0,2))
legend(x = 1.3, y = 320,
       legend = c("NA", "male","female"),
       fill =  rev(sexCol))

##-- NGS samples per sex
statusCol <- c(hcl.colors(3),"gray60")
plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
plot(transects_west, col = "black", add = T, lwd = 0.5)
plot(st_geometry(ngs),
     add = T, pch = 21, cex = 0.5,
     bg = ifelse(ngs$Status == "alpha", statusCol[1],
                 ifelse(ngs$Status == "pup",statusCol[2],
                        ifelse(ngs$Status == "other",statusCol[3],
                               statusCol[4]))),
     col = "black")
legend(  x = 502000, y = 5110000, bty = "n",
         legend = c("alpha", "pup", "other", "NA"),
         title = "NGS samples",
         pch = 21, pt.bg = statusCol)

barplot( cbind(table(status,useNA = "always")),
         col = statusCol,
         beside = F, width = 1, xlim = c(0,2))
legend(x = 1.3, y =320,
       legend = c("NA", "other", "pup", "alpha"),
       fill =  rev(statusCol))


graphics.off()



## -----------------------------------------------------------------------------
## ------ IV. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  betaHab ~ dnorm(0.0,0.01)
  
  habIntensity[1:n.habWindows] <- exp(betaHab*hab.covs[1:n.habWindows])
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
  
  for(i in 1:n.individuals){ 
    z[i] ~ dbern(psi)
  }#i 								
  
  
  ##---- DETECTION PROCESS 
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,6)
  
  betaDet ~ dnorm(0.0,0.01)
  
  logit(p0Traps[1:n.detectors]) <- logit(p0) + betaDet*det.covs[1:n.detectors]
  
  for(i in 1:n.individuals){
    y[i,1:n.maxDets] ~ dbinomLocal_normal(
      size = size[1:n.detectors],
      p0Traps = p0Traps[1:n.detectors],
      sigma = sigma,
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
                      n.habWindows = habitat$n.HabWindows,
                      n.detectors = detectors$n.detectors, 
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      n.maxDets = dim(yCombined.aug)[2],
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

nimData <- list( y = yCombined.aug,
                 z = c(rep(1,n.detected),
                       rep(NA,dim(yCombined.aug)[1]-n.detected)),
                 lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords,
                 size = rep(1,detectors$n.detectors),
                 detCoords = detectors$scaledCoords,
                 hab.covs = c(scale(st_drop_geometry(habitat$grid[ ,c("elev")]))),
                 det.covs = c(scale(st_drop_geometry(detectors$grid[ ,c("transect_L")]))),
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects$habitatGrid)

nimParams <- c("N", "p0", "betaDet", "betaHab", "sigma", "psi", "z", "s")



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
  
  z.init <- rbinom(n = nimConstants$n.individuals, 1, prob = 0.5)
  z.init[!is.na(nimData$z)] <- NA
  
  nimInits <- list( "s" = s.init,
                    "z" = z.init,
                    "psi" = 0.6,
                    "betaDet" = 0,
                    "betaHab" = 0,
                    "p0" = c(0.05),
                    "sigma" = c(1))
  
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
      file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))

pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput)
graphics.off()

plot(nimOutput)




## -----------------------------------------------------------------------------