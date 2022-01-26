################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----- D[...] ------- #####
##### ----- p0[...] ------ #####
##### ----- sigma[] ------ #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(rgdal)
library(raster)
library(sp)
library(coda)
library(nimble)
library(nimbleSCR)
#library(spdep)
library(rgeos)
library(maptools)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
#library(ggplot2)
library(ncdf4)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS (IF NECESSARY) ------
#sourceDirectory(file.path(getwd(),"Source"))

source(file.path(getwd(),"Source/SampleTrack.R"))

## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## Model name
modelName = "AW.0.0"
WD = analysisDir

## HABITAT SPECIFICATIONS
habitat = list( country =  c("SWE","NOR"),
                resolution = 20000, 
                buffer = 60000)

## NGS DATA SPECIFICATIONS
dna = list( sex = c("female","male"),
            year = 2019,
            samplingMonths = list(12,1:6)) 

## DETECTORS SPECIFICATIONS
detectors = list( detSubResolution = 2000,
                  detResolution = 10000,
                  detDeadResolution = 15000)

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(WD, modelName))){dir.create(file.path(WD, modelName))}



## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
## Study area grid
studyAreaGrid <- readOGR(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- spTransform(x = studyAreaGrid, CRSobj = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
studyArea <- aggregate(studyAreaGrid)

## NEED COUNTRY LEVEL POLYGON TO CUT OUT SEA AND OTHER NON-HABITAT
countries <- readOGR(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## Plot check
plot(studyAreaGrid, col = "gray80")
plot(countries, col = "gray80", add = T)
plot(studyAreaGrid, col = adjustcolor("forestgreen",0.2), add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 1, ], col = "gray40", add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 2, ], col = "gray60", add = T)
plot(studyAreaGrid[studyAreaGrid$SCR == 3, ], col = "firebrick3", add = T)



## ------   2. SEARCH EFFORT DATA ------
## Load GPS search transects
transects <- readOGR(file.path(dataDir, "/GISData/Transects_Wolfalps20202021/wolfalps_transects_20202021.shp"))

## Plot check
plot(transects, col = "blue", add = T)



## ------   3. NGS DATA ------
# ngs <- read.csv(file.path(dataDir,
#                           "RovbaseData/ROVBASE DOWNLOAD 20191008/dna_wolf.csv"))
# names(ngs) <- c("ngs_ID", "ngs_code", "RovBaseId", "Origin", "Type",
#                 "Date", "North", "East", "Species", "Analysis_status",
#                 "ID", "Sex", "Commune_num", "Commune", "County_num", "County")
# ## Format dates
# ngs$Date <- as.POSIXct(strptime(ngs$Date, "%d.%m.%Y"))
# ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
# ngs$Month <- as.numeric(format(ngs$Date,"%m"))
# 
# ## Subset to ngs samples from 2019
# ngs <- ngs[ngs$Year %in% dna$year, ]
# 
# ## Remove samples without coordinates
# ngs <- ngs[!is.na(ngs$North), ]
# 
# ## Remove samples without ID
# ngs <- ngs[!is.na(ngs$ID), ]
# 
# ## Transform to SpatialPointsDataFrame
# coordinates(ngs) <- cbind(ngs$East, ngs$North)
# proj4string(ngs) <- proj4string(countries)
# 
# ## Remove samples outside Sweden & Norway
# ngs <- ngs[!is.na(over(ngs, countries)[ ,1]), ]
# 
# ## Plot check
# plot(ngs, add = T)



## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. HABITAT ------
## ------     1.1. HABITAT CHARACTERISTICS ------
## ----- Create Habitat polygon (study area + buffer)  
habitat <- rgeos::gBuffer(studyArea, width = 60000)
## Remove un-suitable habitat (e.g. seas)
habitat <- gIntersection(aggregate(countries), habitat, drop_lower_td = TRUE) 

## ----- Create Habitat raster (keep raster cells with > 50% habitat) 
habResolution <- 10000
habRaster <- raster(extent(habitat))                                         
res(habRaster) <- habResolution/10 
habRaster <- fasterize(st_as_sf(habitat), habRaster)
habRaster[is.na(habRaster)]<- 0
habRaster <- aggregate( habRaster, fact = 10, sum, na.rm = TRUE)
habRaster[habRaster < 50] <- 0
habRaster[habRaster >= 50] <- 1

## ----- Create Habitat matrix of cell ID 
habMatrix <- habRaster
## Identify suitable habitat cells
isHab <- which(habRaster[]==1)
## Cell ID starts from the top left corner, increments left to right and up to down
habMatrix[isHab] <- 1:length(isHab)
## Convert to matrix
habMatrix <- as.matrix(habMatrix)

## ----- Obtain xy coordinates of habitat cells   
habCoords <- coordinates(habRaster)[isHab, ]
dimnames(habCoords) <- list(1:length(isHab), c("x","y"))
# habitat.sp <- SpatialPointsDataFrame(data.frame(habCoords[ ,c("x","y")]),
#                                      data = data.frame(habCoords),
#                                      proj4string = CRS(projection(habitat)))
## Retrieve habitat windows boundaries
loHabCoords <- habCoords - 0.5*habResolution
upHabCoords <- habCoords + 0.5*habResolution
n.HabWindows <- dim(loHabCoords)[1] ## == length(isHab)

## ----- Visual plotting to check if everything is right  
plot(habRaster)					
plot(countries, add=T)
plot(habitat, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
plot(studyArea, add = T, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 1))

## ----- Save mask to crop Environmental Covariates
#crs(habRaster) <- crs(studyAreaGrid)
#writeRaster(habRaster, "mask_habitat", format="GTiff", overwrite=TRUE)



## ------     1.2. HABITAT COVARIATES ------

## ------       1.2.1 CORINE LAND COVER ------

## Load data
CLC <- raster(file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/corine32Ncrop_reclassed.tif"))
# habPoly <- rasterToPolygons(habRaster)
# td <- file.path(dataDir,"/GISData/Environmental Layers/CLC 2018/")
# writeOGR(habPoly, td , "hab_shapefile", driver="ESRI Shapefile", overwrite_layer =TRUE)

grid_hab <- readOGR(file.path(dataDir,"/GISData/shape_studyarea_ALPS/Habitat_shapefile.shp"))

## Plot CLC and the Habitat grid
plot(CLC)
plot(grid_hab, add=T)

#Elements for each class in each cell
#extract for each cell the number of each element per class
# THIS WILL TAKE FOREVER
ext = raster::extract(CLC,grid_hab,df=T) 
# write.csv(ext,"lcl.csv")

#df=true convert  large list in dataframe
#% land calculation: table
prop.clc = ext %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = (n / sum(n))*100) %>%    # add new variable
  # calculate the % of land cover for each class in each cell
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, tidyr::nesting(lc_type), 
                  fill = list(pland = 0)) %>%    # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                  # convert to long format
write.csv(prop.clc,"perc_lcl.csv")
# Legend
# 1	= developed
# 2 = agriculture
# 3 = forest
# 4 = herbaceous
# 5 =	bare rock 
# 6 = sea features
# 7 = inland water
# 9 = perpetual snow

## ------   1.2.2 ELEVATION -------
## Load data
DEM <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/E40N20_crop_1km.tif"))

## Plot DEM and the Habitat grid
plot(DEM)
plot(grid_hab, add=T)

# THIS WILL TAKE FOREVER
ext_el <- raster::extract(DEM, grid_hab, df=T)


mean_el <- ext_el %>%
  setNames(c("ID", "elev")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(elev), n = n())     # calculate mean for cell (ID)


write.csv(mean_el,"mean_el.csv")


## ------   1.2.3 TRI (Terrain Ruggedness Index) -------
## Load data
tri <- raster(file.path(dataDir,"/GISData/Environmental Layers/Orography/TRI_1km.tif"))

## Plot DEM and the Habitat grid
plot(tri)
plot(grid_hab, add=T)

# THIS WILL TAKE FOREVER
ext_tri <- raster::extract(tri, grid_hab, df=T)


mean_tri <- ext_tri %>%
  setNames(c("ID", "tri")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(tri), n = n())     # calculate mean for cell (ID)


write.csv(mean_tri,"mean_tri.csv")

## ------   1.2.4 Human Population Density -------
## Load data
HPop <- raster(file.path(dataDir,"/GISData/Environmental Layers/Human Population/Pop_Dens_1km.tif"))

## Plot DEM and the Habitat grid
plot(HPop)
plot(grid_hab, add=T)

# THIS WILL TAKE FOREVER
ext_hpop <- raster::extract(HPop, grid_hab, df=T)


mean_hpop <- ext_hpop %>%
  setNames(c("ID", "hp")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(hp), n = n())     # calculate mean for cell (ID)


write.csv(mean_hpop,"mean_humpop.csv")


## ------   2. DETECTORS ------
## ------     2.1. DETECTORS CHARACTERISTICS ------
detResolution <- 5000
## Create a raster of searched grid cells
detector.r <- raster(extent(transects),
                     resolution = detResolution,
                     crs = projection(transects))

# ## sf way ----
# rsp <- rasterToPolygons(detector.r)
# rsp <- st_as_sf(rsp)
# transects.sf<- st_as_sf(transects) 
# # intersection
# int = st_intersection(transects.sf, rsp)
# # find out about the length of each line segment
# int$len = st_length(int)
# # use a meaningful id (so far consists only of 0s)
# rsp$Id = 1:nrow(rsp)
# # spatial overlay
# join = st_join(rsp, int)
# # use the ID of the polygon for the aggregation
# out = group_by(join, Id) %>%
#   summarize(length = sum(len))
# # find out about polygons without line segments 
# filter(out, is.na(length))
# # you can set the length of the polygons without line intersections to 0 
# # if you want
# mutate(out, length = ifelse(is.na(length), 0, length))
# det <- st_coordinates(out[!is.na(out$length), ])

## Old version ----
transects.p <- SampleTrack( x = transects,
                            d = 50)
detector.r <- rasterize( transects.p,
                         detector.r,
                         field = "LID",
                         fun = "count")
detector.r[] <- detector.r[]*20
isDet <- which(!is.na(detector.r[])) 
detCoords <- cbind.data.frame( "x" = coordinates(detector.r)[isDet,1],
                               "y" = coordinates(detector.r)[isDet,2],
                               "length" = detector.r[isDet])
plot(detector.r, add=T)
plot(transects, col = "blue", add = T)



## ------     2.2. DETECTOR COVARIATES ------

## ------     2.2.1 SNOWFALL ERA-5 ------


f <- file.path(dataDir,"/GISData/Environmental Layers/Snowfall_2020-2021/ERA-5_snowfall_alps_32N.nc")
oct20 <- brick(f, var="Band1", lvar=4)
nov20 <- brick(f, var="Band2", lvar=4)
dec20 <- brick(f, var="Band3", lvar=4)
jan21 <- brick(f, var="Band4", lvar=4)
feb21 <- brick(f, var="Band5", lvar=4)
mar21<- brick(f, var="Band6", lvar=4)
ap21 <- brick(f, var="Band7", lvar=4)

# ----   OCTOBER 2020  ----

plot(oct20)
plot(transects,col = "blue", add=T)

ext_octsf <- raster::extract(oct20, transects, df=T)
mean_octsf <- ext_octsf %>%
  setNames(c("ID", "sf_oct")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_oct), n = n())     # calculate mean for cell (ID)

# ----   NOVEMBER 2020  ----


plot(nov20)
plot(transects,col = "blue", add=T)

ext_novsf <- raster::extract(nov20, transects, df=T)

mean_novsf <- ext_octsf %>%
  setNames(c("ID", "sf_nov")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_nov), n = n())     # calculate mean for cell (ID)

# ----    DECEMBER 2020 -----

plot(dec20)
plot(transects,col = "blue", add=T)

ext_decsf <- raster::extract(dec20, transects, df=T)

mean_decsf <- ext_decsf %>%
  setNames(c("ID", "sf_dec")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_dec), n = n())     # calculate mean for cell (ID)


# ----   JANUARY 2021  ----

plot(jan21)
plot(transects,col = "blue", add=T)

ext_jansf <- raster::extract(jan21, transects, df=T)

mean_jansf <- ext_jansf %>%
  setNames(c("ID", "sf_jan")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_jan), n = n())     # calculate mean for cell (ID)


# ----   FEBRUARY 2021 ----

plot(feb21)
plot(transects,col = "blue", add=T)

ext_febsf <- raster::extract(feb21, transects, df=T)

mean_febsf <- ext_febsf %>%
  setNames(c("ID", "sf_feb")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_feb), n = n())     # calculate mean for cell (ID)


# ----   MARCH 2021 ----

plot(mar21)
plot(transects,col = "blue", add=T)

ext_marsf <- raster::extract(mar21, transects, df=T)

mean_marsf <- ext_marsf %>%
  setNames(c("ID", "sf_mar")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_mar), n = n())     # calculate mean for cell (ID)



# ----   APRIL 2021 ----

plot(apr21)
plot(transects,col = "blue", add=T)

ext_aprsf <- raster::extract(apr21, transects, df=T)

mean_aprsf <- ext_aprsf %>%
  setNames(c("ID", "sf_apr")) %>%        # rename for ease
  group_by(ID) %>%               # group by ID
  summarise(mean = mean(sf_apr), n = n())     # calculate mean for cell (ID)


## ------   3. DETECTION DATA ------
## ------     3.1. DETECTION MATRIX : y[i,j] ------
# y.ar <- MakeY( myData = myData.alive$myData.sp,
#                myDetectors = myDetectors$main.detector.sp,
#                method = "Binomial",
#                myData2 = myData.dead,
#                myDetectors2 = myDetectors.dead$main.detector.sp,
#                returnIdvector = TRUE)

## SIMULATE DATA FOR NOW?

## ------     3.2. INDIVIDUAL COVARIATES ------






## ------     3.3. DATA AUGMENTATION ------
y <- MakeAugmentation( y = y,
                       aug.factor = augFactor,
                       replace.value = 0)

sex <- MakeAugmentation( y = sex,
                         aug.factor = augFactor,
                         replace.value = NA)


## -----------------------------------------------------------------------------
## ------ III. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##-- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] ~ dnorm(0.0,0.01)
  }
  
  habIntensity[1:n.habWindows] <- exp( betaHab[1:n.habCovs] %*% 
                                          hab.covs[1:n.habWindows,1:n.habCovs])
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i, 1:2, 1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##-- DEMOGRAPHIC PROCESS  
  psi ~ dunif(0,1)   
  rho ~ dunif(0,1)
  
  for(i in 1:n.individuals){ 
    sex[i] ~ dbern(rho)
    z[i] ~ dbern(psi)
  }#i 								
  
  
  ##-- DETECTION PROCESS 
  p0[1] ~ dunif(0,1)
  p0[2] ~ dunif(0,1)

  sigma[1] ~ dunif(0,4)
  sigma[2] ~ dunif(0,4)
  
  for(c in 1:n.detCovs){
    betaDet[c] ~ dnorm(0.0,0.01)
  }#c
  
  for(j in 1:n.detectors){
    logit(p0Vec[j]) <- logit(p0) + 
      betaDet[1:n.detCovs] %*% det.covs[j,1:n.detCovs]
  }#j
  
  for(i in 1:n.individuals){
    p0ID[i] <- p0[sex[i]+1]
    sigmaID[i] <- sigma[sex[i]+1]
    y[i,1:n.maxDetectors,t] ~ dbinomLocal_normal(
      size = 1,
      p0 = p0ID[i],
      p0Traps = ,
      sigma = ,
      s = ,
      trapCoords = ,
      localTrapsIndices = ,
      localTrapsNum = ,
      habitatGrid = ,
      indicator = z[i],
      lengthYCombined = )
  }#i
  
  
  ##---------------------------------------------------------------------------------------------										  
  ##-- DERIVED PARAMETERS 
  N <- sum(z[1:M])
})



## ------   2. CONSTANTS ------
nimConstants <- list( n.individuals = dim(y)[1],
                      n.habWindows = nHabCells,
                      n.habCovs = dim(detCovs)[3],
                      n.detectors = nDetectors,  
                      n.detCovs = max(detCounties),
                      y.max = dim(habIDCells.mx)[1],
                      x.max = dim(habIDCells.mx)[2])



## ------   3. INITIAL VALUES ------
simInits <- list( "z" = z.init$z.init.mx,
                  "p0" = array(runif(18,0,0.2), c(nimConstants$n.countries,dim(y.alive)[3])),
                  "sigma" = runif(1,0,4),
                  "betaDens" = runif(1,-0.1,0.1),#[CM]#0,
                  "betaCovs" = runif(dim(detCovs)[3],-0.1,0.1),#[CM]rep(0,dim(detCovs)[3]),
                  "betaResponse" = runif(1,-0.1,0.1),#[CM]#0,
                  "h" = runif(dim(y.alive)[3]-1,0.1,0.3), #  runif(dim(y.alive)[3]-1,0.2,0.4),
                  "w" = runif(dim(y.alive)[3]-1,0.1,0.3)) #,runif(dim(y.alive)[3]-1,0.2,0.4))



## ------   4. DATA ------
nimData <- list( z = z.init$z.reconstruct.mx,   
                 y.alive = y.alive,
                 y.dead = y.dead,#y.dead.binary,
                 lowerHabCoords = lowerHabCoords/1000, 
                 upperHabCoords = upperHabCoords/1000, 
                 detCountries = detCounties,#detCountries,
                 detCovs = detCovs,
                 detResponse = already.detected ,#+ 1,
                 denCounts = denCounts,
                 trials = n.trials,
                 alpha = rep(1,3),
                 detector.xy = detector.xy/1000,
                 sxy = sxy.data,
                 habitatGrid = habIDCells.mx )

## ------   5. PARAMETERS ------
nimParams <- c("N", 
               "dispSigma",
               "omeg1", "gamma", "phi", "h", "w",
               "p0", "sigma", "betaDens", "betaCovs","betaResponse",
               "z", "sxy")


## ------   6. SAVE NIMBLE INPUT ------
save(nimData,
     nimConstants,
     nimParams,
     modelCode,
     nimInits,
     file = file.path(myVars$WD, myVars$modelName,
                      paste(myVars$modelName, c, ".RData", sep = "")))


apply(nimData$z, 2, function(x) sum(x==2,na.rm = T))



## ------   7. FIT NIMBLE MODEL ------
load(file.path(myVars$WD, myVars$modelName,
               paste(myVars$modelName, c, ".RData", sep = "")))

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 

conf <- configureMCMC(model, monitors = nimParams, thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc



## ------   8. RUN NIMBLE MCMC IN SUCCESSIVE BITES ------
## SET NUMBER OF BITES AND NUMBER OF ITERATIONS PER BITE
bite.size <- 100 
bite.number <- 10

## LOOP OVER NUMBER OF BITES
for(nb in 1:bite.number){
  print(nb)
  if(nb == 1){
    ## run initial MCMC
    MCMCRuntime <- system.time(Cmcmc$run(bite.size))
  } else {      
    ## run subsequent MCMCs
    MCMCRuntime <- system.time(Cmcmc$run(bite.size, reset = FALSE))
  }
  
  ## STORE BITE OUTPUT IN A MATRIX
  mcmcSamples <- as.matrix(Cmcmc$mvSamples)
  CumulRunTime <- proc.time() - ptm
  
  ## EXPORT NIMBLE OUTPUT 
  outname <- file.path(path.OUT, paste("NimbleBite", nb, "_FOR", set, sep = ""))
  save(CumulRuntime, MCMCRuntime, mcmcSamples, file = outname)
  
  ## FREE UP MEMORY SPACE 
  rm("mcmcSamples") 
  Cmcmc$mvSamples$resize(0) ## reduce the internal mvSamples object to 0 rows,
  gc() ## run R's garbage collector
}#nb




## -----------------------------------------------------------------------------
## ------ IV. PROCESS RESULTS ------
## -----------------------------------------------------------------------------
