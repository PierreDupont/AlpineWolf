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

sourceCpp("C:/My_documents/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/cpp/GetDensity_PD.cpp")


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "WesternAlps.4"
thisDir <- file.path(analysisDir, modelName)

## HABITAT SPECIFICATIONS
habitat = list( country =  c("SWE","NOR"),
                resolution = 5000, 
                buffer = 20000)

## DETECTORS SPECIFICATIONS
detectors = list( detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

## NGS DATA SPECIFICATIONS
data = list( sex = c("female","male"),
             status = c("alpha","pup","other"),
             aug.factor = 8) 

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




## ------   2. HABITAT ------
## ------     2.1. HABITAT CHARACTERISTICS ------
##---- Create Habitat polygon (study area + buffer)  
## Subset the study area to western Alps 
# buffData <- st_union(st_buffer(ngs, dist = 20000))

## Buffer the detector grid to delineate the habitat polygon
habitat$polygon <- st_buffer(st_union(detectors$grid), habitat$buffer)


## Remove un-suitable habitat (e.g. seas)
habitat$polygon <- st_intersection( st_union(countries),
                                    habitat$polygon,
                                    drop_lower_td = TRUE)
#plot(habitat$polygon, add=T, col = adjustcolor("gray20",0.2))

##---- Create Habitat raster (keep raster cells with > 50% habitat) 
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
habitat$lowerCoords <- habitat$coords - 0.5*habitat$resolution
habitat$upperCoords <- habitat$coords + 0.5*habitat$resolution
habitat$n.HabWindows <- dim(habitat$lowerCoords)[1] ## == length(isHab)

##---- Visual plotting to check if everything is right  
plot(habitat$raster)					
plot(countries["CNTR_CODE"], col = NA, add = T)
plot(habitat$polygon, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
#plot(studyArea_west, add = T, col = "gray40")
plot(transects,add=T, col="red")
plot(detectors$grid,add=T)





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


# INITIALIZE EMPTY ARRAYS FOR THE DETECTIONS AND DETECTOR INDICES
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


## Plot check
{
  i <- sample(n.detected,1)
  plot(habitat$polygon)
  plot(st_geometry(detectors$grid),add=T)
  plot(detectors$grid[detIndices[i,1:detNums[i]], ],add=T,col = "green")
  plot(ngs[ngs$Genotype.ID == rownames(detMat)[i], ], add = T, col = "red", pch = 3)
}

## ------     4.2. DATA AUGMENTATION ------
yCombined <- MakeAugmentation( y = yCombined,
                               aug.factor = data$aug.factor,
                               replace.value = 0)


## -----------------------------------------------------------------------------
## ------ IV. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
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
  
  for(i in 1:n.individuals){
    y[i,1:n.maxDets] ~ dbinomLocal_normal(
      size = size[1:n.detectors],
      p0 = p0,
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
nimConstants <- list( n.individuals = dim(yCombined)[1],
                      n.habWindows = habitat$n.HabWindows,
                      n.detectors = detectors$n.detectors, 
                      n.localIndicesMax = localObjects$numLocalIndicesMax,
                      n.maxDets = dim(yCombined)[2],
                      y.max = dim(habitat$matrix)[1],
                      x.max = dim(habitat$matrix)[2])

nimData <- list( y = yCombined,
                 z = c(rep(1,n.detected),
                       rep(NA,dim(yCombined)[1]-n.detected)),
                 lowerHabCoords = habitat$loScaledCoords, 
                 upperHabCoords = habitat$upScaledCoords,
                 habIntensity = rep(1,habitat$n.HabWindows),
                 size = rep(1,detectors$n.detectors),
                 detCoords = detectors$scaledCoords,
                 localTrapsIndices = localObjects$localIndices,
                 localTrapsNum = localObjects$numLocalIndices,
                 habitatGrid2 = habitat$matrix,
                 habitatGrid = localObjects$habitatGrid)

nimParams <- c("N", "p0", "sigma", "psi", "z", "s")


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
for(c in 2:3){
  print(system.time(
    runMCMCbites( mcmc = Cmcmc,
                  bite.size = 500,
                  bite.number = 10,
                  path = file.path(thisDir, paste0("output/chain",c))) 
  ))
}#c

##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 0,
                               param.omit = c("s","z"))



##---- Save MCMC output & plot
save( nimOutput,
      file = file.path(thisDir, "mcmc_output.RData"))

pdf(file = file.path(thisDir, "traceplots.pdf"))
plot(nimOutput)
graphics.off()


## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   2. DENSITY ------
## ------     2.1. PREPARE INPUT -----
##-- Load MCMC samples
nimOutput_full <- collectMCMCbites( path = file.path(thisDir, "output"),
                                    burnin = 0)
myResults <- ProcessCodaOutput(nimOutput_full)


## Extract the habitat raster 
habitat.r <- habitat$raster
## Set cells outside the habitat to NA
habitat.r[habitat.r == 0] <- NA
# ## Rasterize to get admin units IDs
# admin.r <- rasterize(Cleanpoly, habitat.r)

## Create a matrix of raster cellsIDs
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = dim(habitat.r)[1],
                      ncol = dim(habitat.r)[2],
                      byrow = TRUE)

## Create a matrix of admin units binary indicators
## (rows == regions ; columns == habitat raster cells)
hab.rgmx <- do.call(rbind, lapply(1, function(x)habitat.r[] == 1))
hab.rgmx[is.na(hab.rgmx)] <- 0
row.names(hab.rgmx) <- "hab"


## ------     2.2. CALCULATE AC-BASED DENSITY ------
RD_Density <- GetDensity_PD(
  sx = myResults$sims.list$s[ , ,1],
  sy = myResults$sims.list$s[ , ,2],
  z = myResults$sims.list$z,
  IDmx = habitat.id,
  aliveStates = 1,
  #ncell = ncell(habitat.r),
  returnPosteriorCells = T,
  regionID = hab.rgmx)

## SAVE AC-based DENSITY FILES
save(RD_Density, file = file.path(thisDir, paste0(modelName,"_density.RData")))




## ------     2.3. DENSITY RASTERS ------
meanDensity.R <- habitat.r
meanDensity.R[] <- RD_Density$MeanCell
plot(meanDensity.R)
meanDensity.R[is.na(habitat.r[])] <- NA
writeRaster( x = meanDensity.R, overwrite = TRUE,
             filename = file.path(thisDir, paste0(modelName, "_density.tif")))



## ------     2.5. DENSITY PLOTS ------
pdf(file = file.path(myVars$WD, myVars$modelName, "PLOTS",
                     paste(myVars$modelName,"_Density.pdf", sep ="")),
    width = 20, height = 15)



## ------       2.5.1. DENSITY MAPS ------
##-- Set color scale
maxDens <- max(RD_Density$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)

##-- Plot overall density raster
meanDensity.R <- habitat.r
meanDensity.R[ ] <- RD_Density$MeanCell
meanDensity.R[is.na(habitat.r[])] <- NA

##-- Overall density plot
par(mfrow = c(1,2), mar = c(4,4,0,0))
plot( habitat$polygon, col = "gray80", border = grey(0.3))
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( habitat$polygon, add = T, border = grey(0.3))
plot(st_geometry(countries), add = T)

mtext( text = paste( "N = ", round(RD_Density$summary["Total",1],1),
                     " [", round(RD_Density$summary["Total",4],1), " ; ",
                     round(RD_Density$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


## ------       2.5.2. PREDICTED DENSITY PLOT ------
##-- Calculate relative density
relativeDens.r <- meanDensity.R
relativeDens.r[] <- meanDensity.R[]/sum(meanDensity.R[], na.rm = T)

##-- Calculate relative point-process intensity
intensity.r <- habitat.r
intensity.r[intensity.r[] > 0] <- exp(
  myResults$mean$betaBNF * nimData$REGION_BNF +
    myResults$mean$betaBNFp * nimData$REGION_BNFp +
    myResults$mean$betaNEU * nimData$REGION_NEU +
    myResults$mean$betaSUMp * nimData$REGION_SUMp +
    myResults$mean$betaR * nimData$R +
    myResults$mean$betaDEM * nimData$DEM +
    myResults$mean$betaDEMR * nimData$DEM * nimData$R)

relativeInt.r <- intensity.r
relativeInt.r[] <- relativeInt.r[]/sum(relativeInt.r[])

##-- Set color scale
maxDens <- max(c(relativeInt.r[], relativeDens.r[]), na.rm = T)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
col <- colFunc(100)

par(mfrow = c(1,2), mar = c(1,1,1,1))
plot( Habitat.list$habitat.poly, border = grey(0.3))
mtext( text = "Realized density", side = 3, line = -4, font = 2)
plot( relativeDens.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( Habitat.list$habitat.poly, add = T)

plot( Habitat.list$habitat.poly, border = grey(0.3))
mtext( text = "Predicted density", side = 3, line = -4, font = 2)
plot( relativeInt.r, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( Habitat.list$habitat.poly, add = T)

## Add density legend
plot( relativeInt.r, breaks = cuts, col = col,
      legend.width = 1, legend.only = TRUE,
      axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 5),
                        labels = round(seq(0,maxDens,length.out = 5)*100, digits = 2), 
                        cex.axis = 0.6),
      legend.args = list( text = 'Relative density (% of pop.km2)', line = -2,
                          side = 4, font = 2, cex = 0.8))

dev.off()



## ------   3. DETECTION ------
pdf(file = file.path(myVars$WD, myVars$modelName, "PLOTS",
                     paste(myVars$modelName,"_Detection.pdf", sep ="")),
    width = 20, height = 15)

par(mfrow = c(1,2))

## ------     3.1. DETECTION MAP ------
##-- DETECTABILITY MAP (p0) 
Detectors.list$main.detector.sp$p0 <- ilogit(logit(myResults$mean$p0[nimData$trapIntIndex]) +
                                               myResults$mean$betaTrap * nimData$DetsCovs)

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
                          lapply(1:length(myResults$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(myResults$sims.list$p0[x,1]) +
                                            myResults$sims.list$betaTrap[x] * detCovs)
                                 }))
p0.posterior2 <- do.call( rbind,
                          lapply(1:length(myResults$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(myResults$sims.list$p0[x,2]) +
                                            myResults$sims.list$betaTrap[x] * detCovs)
                                 }))
p0.posterior3 <- do.call( rbind,
                          lapply(1:length(myResults$sims.list$betaTrap),
                                 function(x){
                                   ilogit(logit(myResults$sims.list$p0[x,3]) +
                                            myResults$sims.list$betaTrap[x] * detCovs)
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
params.simple <- names(myResults$mean)[2:4]
params.means <- do.call(c, lapply(params.simple, function(x)myResults$mean[[x]][1]))
params.sd <- do.call(c, lapply(params.simple, function(x)myResults$sd[[x]][1]))
params.q2.5 <- do.call(c, lapply(params.simple, function(x)myResults$q2.5[[x]][1]))
params.q97.5 <- do.call(c, lapply(params.simple, function(x)myResults$q97.5[[x]][1]))
params.Rhat <- do.call(c, lapply(params.simple, function(x)myResults$Rhat[[x]][1]))
params.n.eff <- do.call(c, lapply(params.simple, function(x)myResults$n.eff[[x]][1]))

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

