## -------------------------------------------------------------------------- ##
## ----------------------- ALPINE WOLF SCR ---------------------------------- ##
## ----------------- s[CLC + HumPop + Zone] --------------------------------- ##
## ---------------- z[psi,rho,theta] ---------------------------------------- ##
## ---------------- y[p0(transects + zone + ),sigma] ------------------------ ##
## -------------------------------------------------------------------------- ##
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(terra)
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
## MODEL NAME 
modelName = "AlpineWolf.5.6_RJMCMC_2324"
thisDir <- file.path(analysisDir, modelName)

## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 30000)

## DETECTORS SPECIFICATIONS
detectors = list( resolution = 5000,
                  detSubResolution = 1000,
                  samplingMonths = c(10:12,1:4))

## NGS DATA SPECIFICATIONS
data = list( sex = c("F","M"),
             status = c("alpha","pup","other"),
             aug.factor = 9) 
# maybe reduce a bit next time 

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}



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

##--- SCR grid
SCRGrid <- read_sf(file.path(dataDir,"GISData/SECR_presence_layer_2020_2021/SECR_presence_layer_2021.shp"))
SCRGrid <- st_transform(x = SCRGrid, crs = st_crs(countries))
SCRGrid <- SCRGrid[SCRGrid$`Pres_20-21` %in% 1, ]

##--- Regions 
regions <- read_sf(file.path(dataDir,"GISData/Output_layout/Alpine_Regions.shp"))
regions <- st_transform(x = regions, crs = st_crs(countries))
regions$ID <- as.numeric(as.factor(regions$DEN_UTS))
# plot(regions)


presence <- read_sf(file.path(dataDir,"GISData/presence_2324/Grid_Alpi_presenzaRAlpine_ITA_2023_2024_LWAEU.shp"))
presence <- st_transform(x = presence, crs = st_crs(countries))
presence <- presence[presence$`2023_2024` %in% 1, ]
presence <- st_intersection(presence, regions)

##--- Alps
alps <- read_sf(file.path(dataDir,"GISData/Output_layout/Italian_Alps.shp"))
alps <- st_transform(x = alps, crs = st_crs(countries))
# plot(alps)


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


## ------   2. SEARCH EFFORT DATA ------
##---- Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20232024/Transetti_def.shp"))

##---- Convert dates
# transects$Date <- parse_date_time(transects$date, orders = c('ymd'))
# transects$Year <- as.numeric(format(transects$Date,"%Y"))
# transects$Month <- as.numeric(format(transects$Date,"%m"))

##---- Plot check
# plot(transects$geometry, col = "red", add = T)

table(st_geometry_type(transects))

# Fix function that replaces invalid lines instead of removing them
# fix_multilinestring <- function(geom) {
#   lines <- st_cast(geom, "LINESTRING", warn = FALSE)
#   # Ensure it's an sfc collection
#   if (!inherits(lines, "sfc")) {
#     lines <- st_sfc(lines, crs = st_crs(geom))
#   }
#   # Fix each LINESTRING
#   fixed_lines <- lapply(lines, function(line) {
#     coords <- st_coordinates(line)[, c("X", "Y")]
#     if (nrow(coords) < 2) {
#       # Duplicate the point to make a minimal valid line
#       coords <- rbind(coords, coords)
#     }
#     st_linestring(coords)
#   })
#   # Build MULTILINESTRING
#   return(st_multilinestring(fixed_lines))
# }
# 
# # Apply to all geometries
# geom_list_fixed <- lapply(st_geometry(transects), fix_multilinestring)
# # Wrap into proper sfc
# fixed_sfc <- st_sfc(geom_list_fixed, crs = st_crs(transects))
# # Combine with original data
# transects_fixed <- st_sf(st_drop_geometry(transects), geometry = fixed_sfc)
# # Track which ones were modified (e.g., had <2-point lines)
# was_fixed <- st_as_text(st_geometry(transects)) != st_as_text(st_geometry(transects_fixed))
# 
# # Add column marking modified features
# transects_fixed$geometry_fixed <- was_fixed
# table(transects_fixed$geometry_fixed)
# 
# 
# fix_multilinestring_safe <- function(geom) {
#   # Force geometry to MULTILINESTRING
#   geom <- st_cast(geom, "MULTILINESTRING", warn = FALSE)
#   # Get coordinates
#   coords_df <- st_coordinates(geom)
#   # If no coordinates, return empty
#   if (nrow(coords_df) == 0) {
#     return(st_geometrycollection())
#   }
#   # Add L1 column if missing
#   if (!"L1" %in% colnames(coords_df)) {
#     coords_df <- cbind(coords_df, L1 = 1)
#   }
#   # Split by sub-line
#   line_list <- split(coords_df[, c("X", "Y")], coords_df[, "L1"])
#   # Fix and ensure all are 2-column matrices
#   fixed_lines <- lapply(line_list, function(mat) {
#     mat <- as.matrix(mat)
#     # Ensure matrix has two columns
#     if (ncol(mat) == 1) {
#       mat <- matrix(mat, ncol = 2)
#     }
#     # Duplicate row if fewer than 2 points
#     if (nrow(mat) < 2) {
#       mat <- rbind(mat, mat)
#     }
#     return(mat)
#   })
#   # Construct MULTILINESTRING safely
#   return(st_multilinestring(fixed_lines))
# }
# 
# fixed_geom_list <- lapply(st_geometry(transects), fix_multilinestring_safe)
# fixed_sfc <- st_sfc(fixed_geom_list, crs = st_crs(transects))
# transects_fixed <- st_sf(st_drop_geometry(transects), geometry = fixed_sfc)
# 
# transects_fixed$geometry_fixed <- st_as_text(st_geometry(transects)) != st_as_text(st_geometry(transects_fixed))
# 
# 
# plot(st_geometry(transects_fixed))
# # Remove/replace non-ASCII characters in all character columns
# transects_fixed_clean <- transects_fixed %>%
#   mutate(across(where(is.character), ~ iconv(., from = "UTF-8", to = "ASCII//TRANSLIT")))
# st_write(transects_fixed_clean, "transects_fixed.shp", delete_layer = TRUE)


##---- All collected NGS samples
# allSamples <- read_sf(file.path(dataDir,"GISData/scats/merged_scats.shp"))
# plot(st_geometry(studyArea))
# plot(st_geometry(countries), col = "gray80",add=T)
# plot(allSamples, add = T,pch=3)


##----CREATE OPERATORS INDEX FOR TRANSECTS -------
# Define the set provinces with less experience
l_provinces <- c("PV", "TV", "AT", "BZ", "TN", "AO", "NO", "BI", "VC", "UD", 
                 "PD", "TS", "GE", "SP", "SV", "SO", "CO", "IM", "BS", "MI", "LC")

# Assign the area_class based on whether the province is in the list
transects$area_class <- ifelse(transects$provinc %in% l_provinces, 1, 2)

institution_index <- c(
  "Aree Protette" = 5,
  "Tecnici" = 5,
  "Forestali" = 4,
  "PolProv" = 4,
  "Universita" = 4,
  "Soccorso" = 4,
  "TecniciCA" = 3,
  "Regione" = 3,
  "Altre Istituzioni" = 2,
  "Volontari" = 1
)

transects$category <- trimws(transects$category)

transects$inst_index <- institution_index[transects$category]
transects$oper_exp <- transects$inst_index * transects$area_class

filtered_transects <- transects %>%
  filter(!REGIONE %in% c("LIG"))

table(filtered_transects$REGIONE)

##---- Genotyped NGS samples -----
ngs <- read.csv(file.path(dataDir,"DNA/ngs_23_24_BL_status.csv"))

dim(ngs)

##---- Filter out scats with no genotypes
ngs <- ngs[!ngs$Genotype.ID == "", ]
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
ngs <- ngs[ngs$Year > 2022, ]
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
ngs$Status[ngs$Status == "pup?"] <- "pup"
ngs$Status[ngs$Status == "alpha?"] <- "alpha"
ngs$Status[ngs$Status == "other?"] <- "other"
ngs$Status[ngs$Status == "unrelated"] <- "other"


table(ngs$Status, useNA = "always")
apply(table(ngs$Genotype.ID, ngs$Status, useNA = "always"),
      2,
      function(x)sum(x)/sum(x>0))

##---- mean number of individual detections by sex and status
apply(table(ngs$Genotype.ID, ngs$Status, ngs$Sex, useNA = "always"),
      c(2,3),
      function(x)sum(x)/sum(x>0))

##---- Turn ngs into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$Coordinate.WGS84.UTM.East, ngs$Coordinate.WGS84.UTM.North)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

##---- Plot check
plot(st_geometry(studyArea))
plot(st_geometry(countries), col = "gray80",add=T)
plot(st_geometry(filtered_transects), col = "red", add = T)
plot(st_geometry(ngs), pch = 3)

# ##---- Calculate distance to the closest transect
# ## For later: can be quite long with many points and transects
# dist <- st_distance(ngs, transects)
# ngs$dist_to_transect <- apply(dist,1,min)
# range(ngs$dist_to_transect)

# ##---- Fix one problematic wolf (detected in 2 distant locations at the same time)
# ngs$Genotype.ID[ngs$Genotype.ID=="WBS-M002"][5:7] <- "WBS-M002.2"
# ngs[ngs$Genotype.ID=="WBS-M002", ]

# ##---- Identify individuals with detections more than 50km apart
# plot(studyArea)
# IDs <- unique(ngs$Genotype.ID)
# resList <- list()
# for( i in 109:length(IDs)){
#   print(i)
#   tmp <- ngs[ngs$Genotype.ID == IDs[i], ]
#   D <- as.numeric(st_distance(tmp))
#   if(any(D > 50000)){
#     resList[[i]] <- tmp
#     print(IDs[i])
#     plot(tmp,pch=19,add=T,col= "red")
#     plot(tmp,pch=19,add=T,type="l",col="red")
#   }
# }
# res <- do.call(rbind,resList)
# write.csv(res,file = "long_distance_detections.csv")


## ------   4. SPATIAL COVARIATES ------
##---- Read snow-fall rasters
SNOW <- rast(file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/snowfall2324.tif"))


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
# plot(CLC)


##---- Load elevation data
DEM_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))
# plot(DEM_raw)


##---- Load terrain ruggedness data
TRI_raw <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))
# plot(TRI_raw)


##---- Load human population density data
POP_raw  <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))
# plot(POP_raw)


##---- Load roads data
# roads <- read_sf(file.path(dataDir,"/GISData/Environmental Layers/Road Density/road_density_tot_clipped.shp"))
# table(roads$highway) # see https://wiki.openstreetmap.org/wiki/Key:highway for road classes description
# mainRoads <- subset(roads, highway %in% c("motorway", "trunk", "primary"))
# rm(roads)


##---- Create a line to delineate the "western Alps"
# cutline <- st_linestring(rbind(c(700000,4800000), c(500000,5230000)))
# cutline <- st_sfc(cutline)
# st_crs(cutline) <- st_crs(studyArea)
# cutline <- st_buffer(cutline, dist = 0.0001)
# plot(studyArea)
# plot(cutline, add = T, border = "blue", lwd = 2)

##---- Create east and western Alps polygons
# studyArea_west <- st_sfc(st_cast(st_difference(st_buffer(studyArea,50000), cutline),"POLYGON")[[1]])
# st_crs(studyArea_west) <- st_crs(studyArea)
# studyArea_east <- st_difference(st_buffer(studyArea,50000),studyArea_west)


##---- Load historical presence data
# packPres_files <- list.files(path = file.path(dataDir, "/GISData/Branchi_SECR"),
#                              pattern = ".shp")
# packPres_raw <- lapply(packPres_files, function(x){
#   tmp <- read_sf(file.path(dataDir, "/GISData/Branchi_SECR", x))
#   tmp <- st_transform(x = tmp, crs = st_crs(studyArea))
#   tmp
# })
# plot(studyArea)
# lapply(packPres_raw, function(x)plot(x, add = T))
# pack.r <- raster(file.path(dataDir, "/GISData/Packs_history_grid/packs_history_sum.tif"))


##---- Load IUCN presence data
list.files(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid"))
iucn_2012_1 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/Clip_2012_12_01_Wolves_permanent.shp"))
iucn_2012_1$SPOIS <- 3
iucn_2012_1 <- st_transform(iucn_2012_1, st_crs(studyArea))
# plot(st_geometry(studyArea))
# plot(st_geometry(iucn_2012_1), add = T, col = "lightskyblue4")

iucn_2012_2 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/Clip_2012_12_01_Wolves_sporadic.shp"))
iucn_2012_2$SPOIS <- 1
iucn_2012_2 <- st_transform(iucn_2012_2, st_crs(studyArea))
# plot(st_geometry(iucn_2012_2), add = T, col = "lightskyblue2")

iucn_2018 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/2018_06_06_Wolf_IUCN_Redlist.shp"))
iucn_2018 <- st_transform(iucn_2018, st_crs(studyArea))
# plot(st_geometry(studyArea))
# plot(iucn_2018[ ,"SPOIS"], col = "blue", add = T)
## Code back to numeric
iucn_2018$SPOIS <- ifelse(iucn_2018$SPOIS == "Sporadic", 1, 3)


iucn_2023 <- read_sf(file.path(dataDir,"/GISData/WOLF_IUCN_LCIE Grid/Wolf_2017-2022_LCIE.shp"))
iucn_2023 <- st_transform(iucn_2023, st_crs(studyArea))
# plot(st_geometry(studyArea))
# plot(iucn_2023[ ,"PRESENCE"], col = "red", add = T)
## Code back to numeric
iucn_2023$SPOIS <- ifelse(iucn_2023$PRESENCE == "Sporadic", 1, 3)

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
# plot(detectors$raster)

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
# plot(detectors$raster)

##---- Remove big cities
to.remove <- POP_raw %>%
  raster::aggregate( x = .,
                     fact = detectors$resolution/res(POP_raw),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>5000) %>%
  st_as_sf() %>%
  fasterize(.,detectors$raster)
detectors$raster[to.remove[ ] == 1] <- NA
# plot(detectors$raster)

##---- Keep original detector area to make habitat
detectors$grid.original <- st_as_sf(rasterToPolygons(detectors$raster))
detectors$grid.original$id <- 1:nrow(detectors$grid.original)
st_crs(detectors$grid.original) <- st_crs(studyArea)

##---- Remove Liguria
to.remove <- regions[regions$ID == 3, ] %>%
  fasterize(.,detectors$raster)
detectors$raster[to.remove[ ] == 1] <- NA
# plot(detectors$raster)

##---- Transform into spatial grid
detectors$grid <- st_as_sf(rasterToPolygons(detectors$raster))
detectors$grid$id <- 1:nrow(detectors$grid)
st_crs(detectors$grid) <- st_crs(studyArea)

##---- Extract length and number of transects in each grid cell
intersection <- st_intersection(detectors$grid, filtered_transects) %>%
  mutate(LEN = .$tt_lngt,
         QI = .$oper_exp) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
            transect_N = sum(n_uscit),   ## Get total number of visits in each detector grid cell
            transect_qi = mean(QI))              ## Get mean transects quality index for each detector grid cell

##---- Store in detector grid
detectors$grid <- detectors$grid %>%
  left_join(intersection, by = "id")
detectors$grid$transect_L[is.na(detectors$grid$transect_L)] <- 0
detectors$grid$transect_N[is.na(detectors$grid$transect_N)] <- 1
detectors$grid$transect_qi[is.na(detectors$grid$transect_qi)] <- 0
detectors$grid$transect_L <- scale(detectors$grid$transect_L)
detectors$grid$mean_transect_L <- scale(detectors$grid$transect_L/detectors$grid$transect_N)

# ##----#---- If you want only grid cells that contain transects instead
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
detectors$size <- as.numeric(table(detectors$sub.sp$main.cell.new.id))

##---- Extract total number of detectors
detectors$n.detectors <- nrow(detectors$grid)



## ------     1.2. DETECTOR COVARIATES ------
## ------       1.2.1. SNOWFALL ERA-5 ------
##---- Calculate mean and cumulative snowfall
# Check the number of layers (should be 6 if it's monthly data for monitpring season)
nlyr(SNOW)
snow_poly <- as.polygons(SNOW, dissolve = FALSE, na.rm = TRUE) %>% st_as_sf()


snow_poly$snow.mean <- rowMeans(st_drop_geometry(snow_poly), na.rm = T)
snow_poly$snow.sum <- rowSums(st_drop_geometry(snow_poly), na.rm = T)

##---- Extract snow-fall in each detector grid cell
intersection <- st_intersection(detectors$grid, snow_poly) %>%
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
  # plot(temp)
  # plot(st_geometry(detectors$grid), add = T)
  
  COV[[l]] <- raster::aggregate( x = temp,
                                 fact = detectors$resolution/res(temp),
                                 fun = mean)
  # plot(COV[[l]])
  # plot(st_geometry(detectors$grid), add = T)
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
# plot(DEM)
# plot(st_geometry(detectors$grid), add = T)

##---- Extract scaled elevation values
detectors$grid$elev <- DEM %>%
  raster::extract(y = detectors$sp) %>%
  scale()




## ------       1.2.4. TERRAIN RUGGEDNESS INDEX -------
##---- Aggregate to the detector resolution
TRI <- raster::aggregate( x = TRI_raw ,
                          fact = detectors$resolution/res(TRI_raw ),
                          fun = mean)

# plot(TRI)
# plot(st_geometry(detectors$grid), add = T)

##---- Extract scaled terrain ruggedness
detectors$grid$tri <- TRI %>%
  raster::extract(y = detectors$sp) %>%
  scale()




## ------       1.2.5. HUMAN POPULATION DENSITY -------
##---- Aggregate to the detector resolution
POP <- raster::aggregate( x = POP_raw ,
                          fact = detectors$resolution/res(POP_raw ),
                          fun = mean)

# plot(POP)
# plot(st_geometry(detectors$grid), add = T)

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
# intersection <- st_intersection(detectors$grid, mainRoads) %>%
#   mutate(length = st_length(.))  %>%
#   st_drop_geometry() %>%
#   group_by(id) %>%
#   summarise(mainRoads_L = sum(length))
# 
# ##---- Store scaled road density
# detectors$grid <- detectors$grid %>%
#   left_join(intersection, by = "id")
# detectors$grid$mainRoads_L[is.na(detectors$grid$mainRoads_L)] <- 0
# detectors$grid$mainRoads_L <- scale(detectors$grid$mainRoads_L)


##---- Display all detector covariates
# plot(detectors$grid[ ,3:length(detectors$grid)], max.plot = 15)


## ------       1.2.7. EAST/WEST ------
##---- Calculate distance of all detectors to each zone
# distWest <- st_distance(detectors$grid,studyArea_west)
# distEast <- st_distance(detectors$grid,studyArea_east)
# 
# ##---- Assign detectors to the closest zone
# detectors$grid$zone <- apply(cbind(distWest,distEast),1,
#                              function(x)which.min(x)-1)
# 
# 


## ------   2. HABITAT ------
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (detector grid + buffer)
habitat$polygon <- st_buffer(st_union(detectors$grid.original),
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
# plot(habitat$raster)

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
# plot(habitat$raster)

##---- Remove big cities
to.remove <- POP_raw %>%
  raster::aggregate( x = .,
                     fact = habitat$resolution/res(POP_raw),
                     fun = mean) %>%
  rasterToPolygons(.,fun = function(x)x>5000) %>%
  st_as_sf() %>%
  fasterize(.,habitat$raster)
habitat$raster[to.remove[ ] == 1] <- NA
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
# plot(habitat$raster)
# plot(st_geometry(countries), add = T)
# plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))
# plot(transects,add=T, col="red")
# plot(st_geometry(detectors$grid), add=T)




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
  # plot(temp)
  # plot(st_geometry(habitat$grid), add = T)
  
  COV[[l]] <- raster::aggregate( x = temp,
                                 fact = habitat$resolution/res(temp),
                                 fun = mean)
  
  COV[[l]] <- raster::focal(x = COV[[l]],
                            w = matrix(1,3,3),
                            fun = mean,
                            na.rm = T)
  
  # plot(COV[[l]])
  # plot(st_geometry(habitat$grid), add = T)
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
# plot(DEM)
# plot(st_geometry(habitat$grid), add = T)

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
# plot(TRI)
# plot(st_geometry(habitat$grid), add = T)

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
# plot(POP)
# plot(st_geometry(habitat$grid), add = T)

## Extract scaled human pop density
habitat$grid$pop <- POP %>%
  raster::extract(y = habitat$sp) %>%
  scale()

## Extract log(human pop density + 1)
log_POP <- POP
log_POP[ ] <- log(POP[]+1)

# par(mfrow=c(1,2))
# plot(POP)
# plot(st_geometry(countries),add=T)
# plot(log_POP)
# plot(st_geometry(countries),add=T)

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
# distWest <- st_distance(habitat$grid, studyArea_west)
# distEast <- st_distance(habitat$grid, studyArea_east)
# 
# ##---- Assign detectors to the closest zone
# habitat$grid$zone <- apply(cbind(distWest,distEast),1,
#                            function(x)which.min(x)-1)
# 

## ------       2.2.7. HISTORICAL PRESENCE ------
# habitat$grid$presence <- 0
# 
# for(y in 1:length(packPres_raw)){
#   ## Extract historical pack presence in each habitat grid cell
#   intersection <- st_intersection(habitat$grid, packPres_raw[[y]]) %>%
#     mutate(pres = 1)  %>%
#     st_drop_geometry() %>%
#     group_by(id) %>%
#     summarise(sumpres = sum(pres))
# 
#   tmp <- habitat$grid %>%
#     left_join(intersection, by = "id")
#   tmp$sumpres[is.na(tmp$sumpres)] <- 0
#   habitat$grid$presence <- habitat$grid$presence + tmp$sumpres
#   plot(habitat$grid[,"presence"])
# }
# 
# habitat$grid$presence <- scale(habitat$grid$presence)
# 



## ------       2.2.8. IUCN PRESENCE ------
## Extract LCIE wolf permanent presence in each habitat grid cell (2012)
intersection <- st_intersection(habitat$grid, iucn_2012_1) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(IUCN = sum(iucn*SPOIS)/2.5e+07)

habitat$grid <- habitat$grid %>%
  left_join(intersection, by = "id")
habitat$grid$IUCN[is.na(habitat$grid$IUCN)] <- 0
# plot(habitat$grid[,"IUCN"])

## Extract LCIE wolf sporadic presence in each habitat grid cell (2012)
intersection <- st_intersection(habitat$grid, iucn_2012_2) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/2.5e+07)

tmp <- habitat$grid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitat$grid$IUCN <- habitat$grid$IUCN + tmp$iucn_2
# plot(habitat$grid[,"IUCN"])

## Extract LCIE wolf presence in each habitat grid cell for 2018 data
intersection <- st_intersection(habitat$grid, iucn_2018) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/2.5e+07)

tmp <- habitat$grid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitat$grid$IUCN <- habitat$grid$IUCN + tmp$iucn_2


## Extract LCIE wolf presence in each habitat grid cell for 2023 data
intersection <- st_intersection(habitat$grid, iucn_2023) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/2.5e+07)

tmp <- habitat$grid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitat$grid$IUCN <- habitat$grid$IUCN + tmp$iucn_2
# plot(habitat$grid[,"IUCN"])

habitat$grid$IUCN <- scale(habitat$grid$IUCN)




## ------       2.2.8. PROTECTED AREAS -----
# intersection <- st_intersection(habitat$grid, PA) %>%
#   mutate(pa = st_area(.))  %>%
#   st_drop_geometry() %>%
#   group_by(id) %>%
#   summarise(PA = sum(pa)/2.5e+07)
#
# ## Store scaled road density
# habitat$grid <- habitat$grid %>%
#   left_join(intersection, by = "id")
# habitat$grid$PA[is.na(habitat$grid$PA)] <- 0
# habitat$grid$PA[habitat$grid$PA > 1] <- 1



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
closest <- nn2( coordinates(detectors$sub.sp),
                st_coordinates(ngs),
                k = 1,
                searchtype = "radius",
                radius = 10000)


##---- Assign each detection to a main detector based on minimum distance
ngs$sub.detector <- c(closest$nn.idx)
ngs$detector <- detectors$sub.sp$main.cell.new.id[closest$nn.idx] 

##---- Drop duplicated detections at the same sub-detectors
ngs <- ngs[!duplicated(ngs[ ,c("sub.detector", "Genotype.ID")]), ] 
#ngs <- droplevels(ngs)

##---- Count individual detections per detector
detMat <- as.matrix(table(ngs$Genotype.ID, ngs$detector))

# ##---- Convert to binary detections (might keep later for binomial model)
# detMat[detMat > 0] <- 1

##---- Retrieve the number of detectors at which each individual is detected
detNums <- apply(detMat, 1, function(x) sum(x>0))

##--- Set-up matrices to store individual detection frequencies and indices
n.detected <- dim(detMat)[1]
detIndices <- matrix(-1, n.detected, max(detNums)*2)
ySparse <- matrix(-1, n.detected, max(detNums)*2)

##---- Fill in the matrix
for(i in 1:n.detected){
  ##-- Get where detections occur (detectors index)
  detIndices[i,1:detNums[i]] <- as.numeric(names(which(detMat[i, ] > 0)))
  ##-- Get detection frequencies (always 1 in case of Bernoulli) 
  ySparse[i,1:detNums[i]] <- detMat[i,which(detMat[i, ] > 0)]
}

##---- Combine individual detection number, frequencies and indices
yCombined <- cbind(detNums, ySparse, detIndices)




## ------     4.2. INDIVIDUAL COVARIATES ------
##---- List detected indivdual names
IDs <- dimnames(yCombined)[[1]]

##---- Set-up vectors to store individual covariates
sex <- status <- pack <- rep(NA, length(IDs))
for(i in 1:length(IDs)){
  ##-- Sex 
  temp <- unique(ngs$Sex[ngs$Genotype.ID %in% IDs[i]])
  if(length(temp)>1)warning(print(paste("ID:", IDs[i], " sex:", temp)))
  sex[i] <- ifelse(length(temp>1),temp[1],temp)
  
  ##- Social status
  temp <- unique(ngs$Status[ngs$Genotype.ID %in% IDs[i]])
  if(length(temp)>1){
    warning(print(paste("ID:", IDs[i], " status:", temp)))
    if(any(temp %in% "na")){status[i] <- NA}
    if(any(temp %in% "na")){status[i] <- "other"}
  } else {
    status[i] <- temp
  }
  
  ##-- Pack membership
  temp <- unique(ngs$Pack[ngs$Genotype.ID %in% IDs[i]])
  if(any(temp %in% "dispersal")){
    status[i] <- temp
  }
}

##---- Convert to numerical values
sex[sex == "F"] <- 0
sex[sex == "M"] <- 1
sex[sex == ""] <- NA
sex <- as.numeric(sex)

status[status == "alpha"] <- 1
status[status == "pup"] <- 2
status[status == "other"] <- 3
status[status == "dispersal"] <- 3
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
{
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
  legend( x = 430000, y = 4900000, bty = "n",
          legend = rev(nameNumDets), title = "# of obs. per detector",
          y.intersp = 1.2, x.intersp = 2,
          pch = 21, pt.bg = rev(colPal[nameNumDets+1]), pt.cex = rev((nameNumDets+1)/2))
  
  
  
  
  ## ------   2. DENSITY COVARIATES ------
  ##-- Corine Land Cover 
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  for(l in 1:length(layer.index)){
    temp <- CLC
    temp[!temp[] %in% layer.index[l]] <- 0
    temp[temp[] %in% layer.index[l]] <- 1
    ## Raw data
    plot(st_geometry(habitat$grid))
    mtext( text = paste("Raw",layer.names[layer.index[l]]), side = 1, font = 2)
    plot(temp, add = T)
    plot(st_geometry(countries), add = T)
    plot(st_geometry(habitat$grid), add = T)
    ## Processed data
    plot(st_geometry(habitat$grid))
    plot(habitat$grid[ ,layer.names[layer.index[l]]], add=T)
    mtext( text = paste("Processed",layer.names[layer.index[l]]), side = 1, font = 2)
    plot(st_geometry(countries), add = T)
    mtext( text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  }#c
  
  # ##-- Elevation
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(habitat$grid))
  # mtext(text = "Raw elevation", side = 1, font = 2)
  # plot(DEM_raw, add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(habitat$grid), add = T)
  # plot(st_geometry(habitat$grid))
  # mtext(text = "Processed elevation", side = 1, font = 2)
  # plot(habitat$grid[ ,"elev"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  # 
  # ##-- Terrain Ruggedness Index
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(habitat$grid))
  # mtext(text = "Raw terrain ruggedness",side = 1, font = 2)
  # plot(TRI_raw, add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(habitat$grid), add = T)
  # plot(st_geometry(habitat$grid))
  # mtext(text = "Processed terrain ruggedness",side = 1, font = 2)
  # plot(habitat$grid[ ,"tri"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  # 
  ##-- Human Population Density 
  par(mfrow = c(1,3), mar = c(2,2,2,2))
  plot(st_geometry(habitat$grid))
  mtext(text = "Raw human density", side = 1, font = 2)
  plot(POP_raw, add = T)
  plot(st_geometry(countries), add = T)
  plot(st_geometry(habitat$grid), add = T)
  plot(st_geometry(habitat$grid))
  mtext(text = "Processed human density", side = 1, font = 2)
  plot(habitat$grid[ ,"pop"], add = T)
  plot(st_geometry(countries), add = T)
  plot(st_geometry(habitat$grid))
  mtext(text = "log(human density+1)", side = 1, font = 2)
  plot(habitat$grid[ ,"log_pop"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  # ##-- Road density
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(habitat$grid))
  # mtext(text = "Raw main roads",side = 1, font = 2)
  # plot(mainRoads[,"highway"], add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(habitat$grid), add = T)
  # plot(st_geometry(habitat$grid[ ,"mainRoads_L"]))
  # mtext(text = "Processed main roads", side = 1, font = 2)
  # plot(habitat$grid[ ,"mainRoads_L"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  # ##-- East/West
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(habitat$grid))
  # mtext(text = "East/West zones",side = 1, font = 2)
  # plot(cutline, add = T, border = "red", lwd = 2)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(habitat$grid[ ,"zone"]))
  # mtext(text = "Processed East/West zones", side = 1, font = 2)
  # plot(habitat$grid[ ,"zone"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- Wolf pack presence
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(st_geometry(habitat$grid))
  mtext(text = "Raw pack presence",side = 1, font = 2)
  lapply(packPres_raw, function(x)plot(x, add = T))
  plot(st_geometry(countries), add = T)
  plot(st_geometry(habitat$grid), add = T)
  plot(st_geometry(habitat$grid[ ,"presence"]))
  mtext(text = "Processed pack presence", side = 1, font = 2)
  plot(habitat$grid[ ,"presence"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- IUCN presence
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(st_geometry(habitat$grid))
  mtext(text = "Raw IUCN presence",side = 1, font = 2)
  plot(iucn_2018, add = T, col = adjustcolor("navyblue",0.5))
  plot(iucn_2012_1, add = T, col = adjustcolor("navyblue",0.2))
  plot(iucn_2012_2, add = T, col = adjustcolor("navyblue",0.5))
  plot(st_geometry(countries), add = T)
  plot(st_geometry(habitat$grid), add = T)
  plot(st_geometry(habitat$grid[ ,"IUCN"]))
  mtext(text = "Processed IUCN presence", side = 1, font = 2)
  plot(habitat$grid[ ,"IUCN"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Density covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  
  
  ## ------   3. DETECTOR COVARIATES ------
  ##-- Search transect length
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  plot(transects, col = "firebrick3", add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Search transects (raw)", side = 1, font = 2)
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  plot(detectors$grid[ ,"transect_L"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Search transects length (processed)", side = 1, font = 2)
  mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- Search transect quality
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  plot(transects, col = "firebrick3", add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Search transects (raw)", side = 1, font = 2)
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  plot(detectors$grid[ ,"transect_qi"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Search transects quality (processed)", side = 1, font = 2)
  mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- Snow
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(habitat$polygon)
  plot(SNOW[ ,"snow.sum"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Cumulative snow fall (raw)", side = 1, font = 2)
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  plot(detectors$grid[,"snow_fall"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Cumulative snow fall(processed)", side = 1, font = 2)
  mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- Corine Land Cover
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  for(l in 1:length(layer.index)){
    temp <- CLC
    temp[!temp[] %in% layer.index[l]] <- 0
    temp[temp[] %in% layer.index[l]] <- 1
    ## Raw data
    plot(studyArea, col = adjustcolor("gray60", alpha.f = 0.3))
    mtext( text = paste("Raw",layer.names[layer.index[l]]), side = 1, font = 2)
    plot(temp, add = T)
    plot(st_geometry(countries), add = T)
    plot(st_geometry(detectors$grid), add = T)
    ## Processed data
    plot(st_geometry(detectors$grid))
    plot(detectors$grid[ ,layer.names[layer.index[l]]], add=T)
    mtext(text = paste("Processed",layer.names[layer.index[l]]), side = 1, font = 2)
    plot(st_geometry(countries), add = T)
    mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  }#c
  
  # ##-- Elevation
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.3))
  # mtext(text = "Raw elevation", side = 1, font = 2)
  # plot(DEM_raw, add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(detectors$grid), add = T)
  # plot(st_geometry(detectors$grid))
  # mtext(text = "Processed elevation", side = 1, font = 2)
  # plot(detectors$grid[ ,"elev"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  # 
  # ##-- Terrain Ruggedness Index
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(detectors$grid))
  # mtext(text = "Raw terrain ruggedness",side = 1, font = 2)
  # plot(TRI_raw, add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(detectors$grid), add = T)
  # plot(st_geometry(detectors$grid))
  # mtext(text = "Processed terrain ruggedness",side = 1, font = 2)
  # plot(detectors$grid[ ,"tri"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- Human Population Density 
  par(mfrow = c(1,3), mar = c(2,2,2,2))
  plot(st_geometry(detectors$grid), add = T)
  plot(st_geometry(detectors$grid))
  mtext(text = "Processed human density", side = 1, font = 2)
  plot(detectors$grid[ ,"pop"], add = T)
  plot(st_geometry(countries), add = T)
  plot(st_geometry(detectors$grid))
  mtext(text = "log(human density+1)", side = 1, font = 2)
  plot(detectors$grid[ ,"log_pop"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  # ##-- Road density
  # par(mfrow = c(1,2), mar = c(2,2,2,2))
  # plot(st_geometry(detectors$grid))
  # mtext(text = "Raw main roads",side = 1, font = 2)
  # plot(mainRoads[,"highway"], add = T)
  # plot(st_geometry(countries), add = T)
  # plot(st_geometry(detectors$grid), add = T)
  # plot(st_geometry(detectors$grid[ ,"mainRoads_L"]))
  # mtext(text = "Processed main roads", side = 1, font = 2)
  # plot(detectors$grid[ ,"mainRoads_L"], add = T)
  # plot(st_geometry(countries), add = T)
  # mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  ##-- East/West
  par(mfrow = c(1,2), mar = c(2,2,2,2))
  plot(st_geometry(detectors$grid))
  mtext(text = "East/West zones",side = 1, font = 2)
  plot(cutline, add = T, border = "red", lwd = 2)
  plot(st_geometry(countries), add = T)
  plot(st_geometry(detectors$grid[ ,"zone"]))
  mtext(text = "Processed East/West zones", side = 1, font = 2)
  plot(detectors$grid[ ,"zone"], add = T)
  plot(st_geometry(countries), add = T)
  mtext(text = "Detector covariates", side = 3, line = -2, outer = TRUE, font = 2)
  
  
  ## ------   4. INDIVIDUAL COVARIATES -----
  ##-- individual centroids
  par(mfrow = c(1,1))
  centroids <- st_drop_geometry(ngs) %>%
    group_by(Genotype.ID) %>%
    summarise(X = mean(CoordX),               ## Get total length searched in each detector grid cell
              Y = mean(CoordY)) 
  coordinates(centroids) <- cbind.data.frame(centroids$X,centroids$Y)
  proj4string(centroids) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
  centroids <- st_as_sf(centroids)
  plot(habitat$polygon, col = "gray20", border = F)
  plot(transects, col = "gray80", add = T)
  plot(centroids, add = T, pch = 21, cex = 1,
       bg = adjustcolor("red",0.2),
       col = "red")
  mtext( text = "Individual centroids", side = 3, font = 2)
  
  ##-- NGS samples per sex
  par(mfrow = c(1,2))
  sexCol <- c(hcl.colors(2),"gray60")
  plot(habitat$polygon, col = adjustcolor("gray60", alpha.f = 0.5), border = F)
  plot(transects, col = "black", add = T)
  plot(st_geometry(ngs), add = T, pch = 21, cex = 1,
       bg = ifelse(ngs$Sex == "F", sexCol[1], sexCol[2]),
       col = "black")
  legend( x = 502000, y = 5110000, bty = "n",
          legend = c("Females", "males", "NA"),
          title = "NGS samples",
          pch = 21,
          pt.bg = sexCol)
  
  par(mar = c(10,10,10,10))
  barplot( cbind(table(sex, useNA = "always")),
           col = sexCol,
           ylab = "Number of individuals detected",
           beside = T)
  legend( "topright",
          legend = c("NA", "male","female"),
          fill =  rev(sexCol))
  
  ##-- NGS samples per sex
  par(mfrow = c(1,2), mar = c(2,2,2,2))
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
           beside = T)
  legend("topright",
         legend = c("alpha","pup","other","NA"),
         fill =  statusCol)
  
  graphics.off()
  
}


## -----------------------------------------------------------------------------
## ------ IV. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  psiRJ ~ dunif(0, 1) # inclusion prob
  
  for(c in 1:n.habCovs){
    betaHab.raw[c] ~ dnorm(0.0,0.01)
    zRJ[c] ~ dbern(psiRJ)
    betaHab[c] <- betaHab.raw[c] * zRJ[c]
  }#c
  
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% (betaHab[1:n.habCovs]*zRJ[1:n.habCovs]))
  
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
    betaDet[c] ~ dnorm(0.0,0.01)
  }
  
  for(st in 1:n.states){
    for(ss in 1:2){
      p0[st,ss] ~ dunif(0,0.5)
      sigma[st,ss] ~ dunif(0,12)
      logit(p0Traps[st,ss,1:n.detectors]) <- logit(p0[st,ss]) + 
        det.covs[1:n.detectors,1:n.detCovs] %*% betaDet[1:n.detCovs]
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
nimData <- list( y = yCombined.aug,
                 z = c(rep(1,n.detected),
                       rep(NA,dim(yCombined.aug)[1]-n.detected)),
                 sex = sex.aug,
                 status = status.aug,
                 alpha = matrix(1,3,2),
                 lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords, 
                 hab.covs = cbind.data.frame(
                   "bare rock" = habitat$grid$`bare rock`,
                   "herbaceous" = habitat$grid$`herbaceous`,
                   "forest" = habitat$grid$`forest`,
                   "pop" = habitat$grid$pop,
                   "IUCN" = habitat$grid$`IUCN`),
                 det.covs = cbind.data.frame(
                   "transect_L" = detectors$grid$`transect_L`,
                   "transect_qi" = detectors$grid$`transect_qi`,
                   "snow_fall" = detectors$grid$`snow_fall`,
                   # "zone" = detectors$grid$`zone`,
                   "log_pop" = detectors$grid$`log_pop`),
                 size = detectors$size,
                 detCoords = detectors$scaledCoords,
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects$habitatGrid)



nimConstants <- list( n.individuals = nrow(nimData$y),
                      n.maxDets = ncol(nimData$y),
                      n.habWindows = habitat$n.HabWindows,
                      n.detectors = detectors$n.detectors, 
                      n.habCovs = ncol(nimData$hab.covs),
                      n.detCovs = ncol(nimData$det.covs),
                      n.states = 3,
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

nimParams <- c("N", "p0", "sigma", "psi", "zRJ","psiRJ",
               "betaDet", "betaHab.raw", "betaHab", "theta", "rho",
               "z", "s", "status", "sex")


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
                    "psiRJ" = 1,
                    "zRJ" = rep(1,nimConstants$n.habCovs),
                    "betaHab.raw" = rep(0,nimConstants$n.habCovs),
                    "p0" = cbind(c(0.1,0.1,0.05),
                                 c(0.1,0.1,0.05)),
                    "sigma" = cbind(c(1,1,2),
                                    c(1,1,2)))
  
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
nimModel$calculate()

##---- Compile the nimble model object to C++
CsimModel <- compileNimble(nimModel)
CsimModel$calculate()

##---- Configure and compile the MCMC object 
conf <- configureMCMC( model = nimModel,
                       monitors = nimParams,
                       thin = 1)

##---- Configure reversible jump
configureRJ( conf =  conf,
             targetNodes = 'betaHab.raw',
             indicatorNodes = 'zRJ',
             control = list(mean = 0, scale = .2))

Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
for(c in 1:1){
  print(system.time(
    runMCMCbites( mcmc = Cmcmc,
                  bite.size = 500,
                  bite.number = 10,
                  path = file.path(thisDir, paste0("output/chain",c)))
  ))
}





## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
##---- Collect multiple MCMC bites and chains
nimOutput_noZ <- collectMCMCbites( path = file.path(thisDir, "output"),
                                   burnin = 25,
                                   param.omit = c("s","z","sex","status"))
res <- ProcessCodaOutput(nimOutput_noZ)

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput_noZ)
graphics.off()


##---- Process and save MCMC samples
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 20,
                               param.omit = nimParams[!nimParams %in%  c("s","z","sex","status")])
nimOutput <- as.mcmc.list(lapply(nimOutput, function(x)as.mcmc(x[seq(1,nrow(x),by=10),])))
gc()
res_sxy <- ProcessCodaOutput(nimOutput)


# res_sxy <- ProcessCodaOutput(nimOutput$samples2)
save(res, res_sxy,
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))

## ------   1. CALCULATE AC-BASED DENSITY ------
##---- Load processed MCMC samples
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
load(file.path(thisDir, "input", paste0(modelName, "_1.RData")))

## ----- CALCULATE DENSITIES WITH A 1X1 KM RESOLUTION -----
##---- Disaggregate habitat raster to a 1x1km resolution
habitat.r <- disaggregate(x = habitat$raster, fact = 5)
italia.r <- disaggregate(x = habitat$Italia, fact = 5)

##---- Create a matrix of raster cellsIDs (including cells outside the habitat)
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = nrow(habitat.r),
                      ncol = ncol(habitat.r),
                      byrow = TRUE)

##---- Rescale sxy coords to original projection and reproject to the new raster
dimnames(res_sxy$sims.list$s) <- list(1:dim(res_sxy$sims.list$s)[1],
                                      1:dim(res_sxy$sims.list$s)[2],
                                      c("x","y"))
s.original <- scaleCoordsToHabitatGrid(
  coordsData = res_sxy$sims.list$s,
  coordsHabitatGridCenter = coordinates(habitat$raster),
  scaleToGrid = F)$coordsDataScaled

##---- ... and reproject to the new raster
s.rescaled <- scaleCoordsToHabitatGrid(
  coordsData = s.original,
  coordsHabitatGridCenter = coordinates(habitat.r),
  scaleToGrid = T)$coordsDataScaled


##---- Create a matrix of Italian regions
##---- (rows == regions ; columns == habitat raster cells)
regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = NA)
regions.r[regions.r[]==0] <- NA
regions.r <- regions.r + italia.r - 1
# for(i in 1:length(regions$DEN_UTS)){
#   regions.r[regions.r %in% regions$ID[i]] <- regions$DEN_UTS[i]
# }
regions.unique <- sort(na.omit(unique(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions.unique

# ##---- Create a matrix of Italian regions without Liguria
# regionsnolig <- regions[regions$DEN_UTS != 'Liguria',]
# 
# regionsnolig.r <- fasterize(sf = st_as_sf(regionsnolig),
#                             raster = habitat.r,
#                             field = "ID",
#                             background = NA)
# 
# regionsnolig.r[regionsnolig.r[]==0] <- NA
# regionsnolig.r <- regionsnolig.r + italia.r - 1
# 
# # for(i in 1:length(regionsnolig$DEN_UTS)){
# #   regionsnolig.r[regionsnolig.r %in% regionsnolig$ID[i]] <- regionsnolig$DEN_UTS[i]
# # }
# regionsnolig.unique <- sort(na.omit(unique(regionsnolig.r[])))
# regionsnolig.rgmx <- do.call(rbind, lapply(regionsnolig.unique, function(x){regionsnolig.r[] == x}))
# regionsnolig.rgmx[is.na(regionsnolig.rgmx)] <- 0
# row.names(regionsnolig.rgmx) <- regionsnolig.unique
# 
# 
# ##---- Create a matrix of grid minimum presence
# ##---- (rows == regions ; columns == habitat raster cells)
# # presence <- st_buffer(st_union(presence), dist = 1000)
# presence.r <- fasterize(sf = st_as_sf(presence),
#                         raster = habitat.r,
#                         field = "ID",
#                         background = NA)
# 
# 
# presence.r[presence.r[]==0] <- NA
# presence.r <- presence.r + italia.r - 1
# # for(i in 1:length(presence$DEN_UTS)){
# #   presence.r[presence.r %in% presence$ID[i]] <- presence$DEN_UTS[i]
# # }
# presence.unique <- sort(na.omit(unique(presence.r[])))
# presence.rgmx <- do.call(rbind, lapply(presence.unique, function(x){presence.r[] == x}))
# presence.rgmx[is.na(presence.rgmx)] <- 0
# row.names(presence.rgmx) <- presence.unique
# 
# 
# # THEN WITHOUT LIGURIA
# presence.nolig <- presence[presence$DEN_UTS != 'Liguria',]
# 
# presencenolig.r <- fasterize(sf = st_as_sf(presence.nolig),
#                              raster = habitat.r,
#                              field = "ID",
#                              background = NA)
# 
# presencenolig.r[presencenolig.r[]==0] <- NA
# presencenolig.r <- presencenolig.r + italia.r - 1
# # for(i in 1:length(presence.nolig$DEN_UTS)){
# #   presencenolig.r[presencenolig.r %in% presence.nolig$ID[i]] <- presence.nolig$DEN_UTS[i]
# # }
# presence.nolig.unique <- sort(na.omit(unique(presencenolig.r[])))
# presence.nolig.rgmx <- do.call(rbind, lapply(presence.nolig.unique, function(x){presencenolig.r[] == x}))
# presence.nolig.rgmx[is.na(presence.nolig.rgmx)] <- 0
# row.names(presence.nolig.rgmx) <- presence.nolig.unique

##---presence.r##---- Create a matrix of the Alpine region
##---- (rows == regions ; columns == habitat raster cells)
alps.r <- fasterize( sf = st_as_sf(alps),
                     raster = regions.r,
                     background = NA)
# plot(alps.r)

alps.r[alps.r[]==0] <- NA
alps.r <- alps.r + italia.r - 1
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"


iter <- seq(1,dim(res_sxy$sims.list$z)[1],length.out = 1000)

##---- Calculate density
WA_Density <- GetDensity(
  sx = s.rescaled[iter, ,1],
  sy = s.rescaled[iter, ,2],
  z = res_sxy$sims.list$z[iter, ],
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = F,
  regionID = regions.rgmx)

gc()
##---- Calculate sex and status-specific densities
WA_status <- list()
for(s in 0:1){
  WA_status[[s+1]] <- list()
  for(ss in 1:3){
    thisStatus <- (res_sxy$sims.list$z[iter, ] == 1) &
      (res_sxy$sims.list$sex[iter, ] == s) &
      (res_sxy$sims.list$status[iter, ] == ss)
    
    WA_status[[s+1]][[ss]] <- GetDensity(
      sx = s.rescaled[iter, ,1],
      sy = s.rescaled[iter, ,2],
      z = thisStatus,
      IDmx = habitat.id,
      aliveStates = 1,
      returnPosteriorCells = F,
      regionID = regions.rgmx)
  }#ss
}#s

gc()
##---- Calculate total and regional densities
# WA_regionsNoLig <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res_sxy$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = regionsnolig.rgmx)
# 
# ##---- Calculate sex and status-specific densities
# WA_status_reg_NoLig <- list()
# for(s in 0:1){
#   WA_status_reg_NoLig[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res_sxy$sims.list$z[iter, ] == 1) &
#       (res_sxy$sims.list$sex[iter, ] == s) &
#       (res_sxy$sims.list$status[iter, ] == ss)
#     
#     WA_status_reg_NoLig[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter, ,1],
#       sy = s.rescaled[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = regionsnolig.rgmx)
#   }#ss
# }#s

# gc()
# # ##---- Calculate overall density
# WA_presence <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = presence.rgmx)
# 
# ##---- Calculate sex and status-specific densities
# WA_presence_status <- list()
# for(s in 0:1){
#   WA_presence_status[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res$sims.list$z[iter, ] == 1) &
#       (res$sims.list$sex[iter, ] == s) &
#       (res$sims.list$status[iter, ] == ss)
#     
#     WA_presence_status[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter, ,1],
#       sy = s.rescaled[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = presence.rgmx)
#   }#ss
# }#s
# ## ------   1.1. EXTRACT PRESENCE ONLY & REGION-SPECIFIC DENSITIES ------
# ##---- Calculate total and regional densities
# WA_presence <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res_sxy$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = presence.rgmx)
# 
# ##---- Calculate sex and status-specific densities
# WA_presence_status <- list()
# for(s in 0:1){
#   WA_presence_status[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res$sims.list$z[iter, ] == 1) &
#       (res$sims.list$sex[iter, ] == s) &
#       (res$sims.list$status[iter, ] == ss)
#     
#     WA_presence_status[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter, ,1],
#       sy = s.rescaled[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = presence.rgmx)
#   }#ss
# }#s
# 
# 
# ##---- Calculate Italian density with no Liguria
# # rownames(presence.rgmx)
# # # Create a logical mask: TRUE when KEEP
# # keep <- !(trimws(tolower(rownames(presence.rgmx))) == "liguria")  
# # # Subsample the object while keeping matrix type
# # presence.nol.rgmx <- as.matrix(presence.rgmx[keep, ])
# 
# WA_presence_noLig <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res_sxy$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = presence.nolig.rgmx)
# 
# # save WA_Italy_noLig and remove the object and reload it
# # save WA_status_noLig for each status and sex and save it remove the object and 
# WA_status_noLig <- list()
# for(s in 0:1){
#   WA_status_noLig[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res_sxy$sims.list$z == 1) &
#       (res_sxy$sims.list$sex == s) &
#       (res_sxy$sims.list$status == ss)
#     
#     WA_status_noLig[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter, ,1],
#       sy = s.rescaled[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = presence.nolig.rgmx)
#   }
# }



## ------   1.2. SAVE 1X1 KM DENSITIES ------
save(WA_Density,
     # WA_status,
     # WA_regionsNoLig,
     # WA_status_reg_NoLig,
     # WA_presence,
     # WA_presence_status,
     # WA_presence_noLig,
     # WA_status_noLig,
     file = file.path(thisDir, paste0(modelName, "_densities_1x1.RData")))
## ------   1.3. LOAD 1X1 KM DENSITIES ------
# load(file.path(thisDir, paste0(modelName, "_densities_1x1.RData")))

## ----  CALCULATE DENSITIES WITH A 5X5 KM RESOLUTION ----

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

regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = NA)
regions.r[regions.r[]==0] <- NA
regions.r <- regions.r + habitat$Italia - 1
for(i in 1:length(regions$DEN_UTS)){
  regions.r[regions.r %in% regions$ID[i]] <- regions$DEN_UTS[i]
}
regions.unique <- sort(na.omit(unique(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions.unique

##---- Create a matrix of Italian regions without Liguria
# regionsnolig <- regions[regions$DEN_UTS != 'Liguria',]
# 
# regionsnolig.r <- fasterize(sf = st_as_sf(regionsnolig),
#                             raster = habitat.r,
#                             field = "ID",
#                             background = NA)
# 
# regionsnolig.r[regionsnolig.r[]==0] <- NA
# regionsnolig.r <- regionsnolig.r + habitat$Italia - 1
# 
# # for(i in 1:length(regionsnolig$DEN_UTS)){
# #   regionsnolig.r[regionsnolig.r %in% regionsnolig$ID[i]] <- regionsnolig$DEN_UTS[i]
# # }
# regionsnolig.unique <- sort(na.omit(unique(regionsnolig.r[])))
# regionsnolig.rgmx <- do.call(rbind, lapply(regionsnolig.unique, function(x){regionsnolig.r[] == x}))
# regionsnolig.rgmx[is.na(regionsnolig.rgmx)] <- 0
# row.names(regionsnolig.rgmx) <- regionsnolig.unique
# 
# 
# ##---- Create a matrix of grid minimum presence
# ##---- (rows == regions ; columns == habitat raster cells)
# # presence <- st_buffer(st_union(presence), dist = 1000)
# presence.r <- fasterize(sf = st_as_sf(presence),
#                         raster = habitat.r,
#                         field = "ID",
#                         background = NA)
# 
# presence.r[presence.r[]==0] <- NA
# presence.r <- presence.r + habitat$Italia - 1
# for(i in 1:length(presence$DEN_UTS)){
#   presence.r[presence.r %in% presence$ID[i]] <- presence$DEN_UTS[i]
# }
# presence.unique <- sort(na.omit(unique(presence.r[])))
# presence.rgmx <- do.call(rbind, lapply(presence.unique, function(x){presence.r[] == x}))
# presence.rgmx[is.na(presence.rgmx)] <- 0
# row.names(presence.rgmx) <- presence.unique

# THEN WITHOUT LIGURIA
# presence.nolig <- presence[presence$DEN_UTS != 'Liguria',]
# 
# presencenolig.r <- fasterize(sf = st_as_sf(presence.nolig),
#                              raster = habitat.r,
#                              field = "ID",
#                              background = NA)
# 
# presencenolig.r[presencenolig.r[]==0] <- NA
# presencenolig.r <- presencenolig.r + habitat$Italia - 1
# # for(i in 1:length(presence.nolig$DEN_UTS)){
# #   presencenolig.r[presencenolig.r %in% presence.nolig$ID[i]] <- presence.nolig$DEN_UTS[i]
# # }
# presence.nolig.unique <- sort(na.omit(unique(presencenolig.r[])))
# presence.nolig.rgmx <- do.call(rbind, lapply(presence.nolig.unique, function(x){presencenolig.r[] == x}))
# presence.nolig.rgmx[is.na(presence.nolig.rgmx)] <- 0
# row.names(presence.nolig.rgmx) <- presence.nolig.unique

##---presence.r##---- Create a matrix of the Alpine region
##---- (rows == regions ; columns == habitat raster cells)
alps.r <- fasterize( sf = st_as_sf(alps),
                     raster = regions.r,
                     background = NA)
# plot(alps.r)

alps.r[alps.r[]==0] <- NA
alps.r <- alps.r + habitat$Italia - 1
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"


iter <- seq(1,dim(res_sxy$sims.list$z)[1],length.out = 1500)

##---- Calculate density
WA_Density <- GetDensity(
  sx = res_sxy$sims.list$s[iter, ,1],
  sy = res_sxy$sims.list$s[iter, ,2],
  z = res_sxy$sims.list$z[iter, ],
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = F,
  regionID = regions.rgmx)

gc()
##---- Calculate sex and status-specific densities
WA_status <- list()
for(s in 0:1){
  WA_status[[s+1]] <- list()
  for(ss in 1:3){
    thisStatus <- (res_sxy$sims.list$z[iter, ] == 1) &
      (res_sxy$sims.list$sex[iter, ] == s) &
      (res_sxy$sims.list$status[iter, ] == ss)
    
    WA_status[[s+1]][[ss]] <- GetDensity(
      sx = res_sxy$sims.list$s[iter, ,1],
      sy = res_sxy$sims.list$s[iter, ,2],
      z = thisStatus,
      IDmx = habitat.id,
      aliveStates = 1,
      returnPosteriorCells = F,
      regionID = regions.rgmx)
  }#ss
}#s

gc()
##---- Calculate total and regional densities
# WA_regionsNoLig <- GetDensity(
#   sx = res_sxy$sims.list$s[iter, ,1],
#   sy = res_sxy$sims.list$s[iter, ,2],
#   z = res_sxy$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = regionsnolig.rgmx)
# 
# gc()
# ##---- Calculate sex and status-specific densities
# WA_status_reg_NoLig <- list()
# for(s in 0:1){
#   WA_status_reg_NoLig[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res_sxy$sims.list$z[iter, ] == 1) &
#       (res_sxy$sims.list$sex[iter, ] == s) &
#       (res_sxy$sims.list$status[iter, ] == ss)
#     
#     WA_status_reg_NoLig[[s+1]][[ss]] <- GetDensity(
#       sx = res_sxy$sims.list$s[iter, ,1],
#       sy = res_sxy$sims.list$s[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = regionsnolig.rgmx)
#   }#ss
# }#s
# 
# gc()
# 
# ###---- Calculate overall density
# WA_presence <- GetDensity(
#   sx = res_sxy$sims.list$s[iter, ,1],
#   sy = res_sxy$sims.list$s[iter, ,2],
#   z = res$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = presence.rgmx)
# 
# gc()
# ##---- Calculate sex and status-specific densities
# WA_presence_status <- list()
# for(s in 0:1){
#   WA_presence_status[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res$sims.list$z[iter, ] == 1) &
#       (res$sims.list$sex[iter, ] == s) &
#       (res$sims.list$status[iter, ] == ss)
#     
#     WA_presence_status[[s+1]][[ss]] <- GetDensity(
#       sx = res_sxy$sims.list$s[iter, ,1],
#       sy = res_sxy$sims.list$s[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = presence.rgmx)
#   }#ss
# }#s
# 
# gc()

###---- Calculate overall density
WA_alps <- GetDensity(
  sx = res_sxy$sims.list$s[iter, ,1],
  sy = res_sxy$sims.list$s[iter, ,2],
  z = res_sxy$sims.list$z[iter, ],
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = F,
  regionID = alps.rgmx)

gc()

## ------   1.2. SAVE 5X5 KM DENSITIES ------
save(WA_Density,
     # WA_status,
     # WA_regionsNoLig,
     # WA_status_reg_NoLig,
     # WA_presence,
     # WA_presence_status,
     # WA_presence_noLig,
     # WA_status_noLig,
     WA_alps,
     file = file.path(thisDir, paste0(modelName, "_densities_5x5.RData")))
## ------   1.3. LOAD 5X5 KM DENSITIES ------
load(file.path(thisDir, paste0(modelName, "_densities_5x5.RData")))


##---- Plot overall density raster
par(mfrow = c(1,1), mar = c(4,4,0,0), bg = "white")
meanDensity.R <- regions.r
meanDensity.R[ ] <- WA_Density$MeanCell
meanDensity.R[is.na(regions.r[])] <- NA
plot(meanDensity.R)


##---- Smooth using a moving window 
f.rast <- function(x) ifelse(is.na(x[13]), NA, mean(x,na.rm = T)) 
MovingWindowSize <-  matrix(1,5,5)
meanDensity.R <- focal(meanDensity.R, MovingWindowSize, f.rast)
plot(meanDensity.R)

# writeRaster(x = meanDensity.R, overwrite = TRUE,
#              filename = file.path(thisDir, paste0(modelName, "regions_raster.tif")))



##---- Plot 1
##---- Set color scale
maxDens <- max(meanDensity.R[],na.rm = T)#max(WA_regions$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   


## ------   2. DENSITY ------
pdf(file = file.path(thisDir, paste0(modelName,"_NoLigDet_morebites_Density.pdf")),
    width = 20, height = 15)

## ------     2.1. DENSITY MAP ------
##---- Set color scale
maxDens <- max(WA_Density$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("whitesmoke", "lightskyblue3", "dodgerblue3","mediumpurple4","lightsalmon1"))
col <- colFunc(100)
par(mfrow = c(2,2), mar = c(4,4,0,0))

##---- Plot overall density raster
meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_Density$MeanCell
meanDensity.R[is.na(habitat.r[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = T)
plot( habitat$polygon, add = T, border = grey(0.3))
plot( st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Density$summary["Total",1],1),
                     " [", round(WA_Density$summary["Total",4],1), " ; ",
                     round(WA_Density$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot Italian density raster
ital.R <- regions.r
ital.R[] <- WA_regions$MeanCell
ital.R[is.na(regions.r[])] <- NA

plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( ital.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = T)
plot( habitat$polygon, add = T, border = grey(0.3))
plot( st_geometry(countries), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_regions$summary["Total",1],1),
                     " [", round(WA_regions$summary["Total",4],1), " ; ",
                     round(WA_regions$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)

# Plot Italia senza assenza 
pres.R <- presence.r
pres.R[] <- WA_presence$MeanCell
pres.R[is.na(presence.r[])] <- NA

plot(habitat$polygon, col = "gray80", border = grey(0.3))   # sfondo: poligono Italia
plot(pres.R, add = TRUE, breaks = cuts, col = col,     # raster di densità Italia - Liguria
     axes = FALSE, box = FALSE, bty = "n", legend = FALSE)
plot(habitat$polygon, add = TRUE, border = grey(0.3))       # bordo Italia
plot(st_geometry(countries), add = TRUE, lwd = 2)            # confini paesi


mtext( text = paste( "N = ", round(WA_presence$summary["Total",1],1),
                     " [", round(WA_presence$summary["Total",4],1), " ; ",
                     round(WA_presence$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)

# Plot Italia presence senza Liguria
presnolig.R <- presencenolig.r
presnolig.R[] <- WA_presence_noLig$MeanCell
presnolig.R[is.na(presencenolig.r[])] <- NA

plot(habitat$polygon, col = "gray80", border = grey(0.3))   # sfondo: poligono Italia
plot(presnolig.R, add = TRUE, breaks = cuts, col = col,     # raster di densità Italia - Liguria
     axes = FALSE, box = FALSE, bty = "n", legend = FALSE)
plot(habitat$polygon, add = TRUE, border = grey(0.3))       # bordo Italia
plot(st_geometry(countries), add = TRUE, lwd = 2)            # confini paesi


mtext( text = paste( "N = ", round(WA_presence_noLig$summary["Total",1],1),
                     " [", round(WA_presence_noLig$summary["Total",4],1), " ; ",
                     round(WA_presence_noLig$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)

graphics.off()

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
pdf(file = file.path(thisDir, paste0(modelName,"regions_Abundance.pdf")),
    width = 20, height = 15)

par(mfrow = c(1,1))
plot.new()
grid.table(round(WA_Density$summary,1))
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

graphics.off()
