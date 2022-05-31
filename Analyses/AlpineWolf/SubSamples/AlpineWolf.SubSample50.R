###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
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
library(gridExtra)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.50"
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
             aug.factor = 6) 

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
SCRGrid <- SCRGrid[SCRGrid$Pres_20.21 == 1, ]

##--- Regions 
regions <- read_sf(file.path(dataDir,"GISData/Output_layout/Alpine_Regions.shp"))
regions <- st_transform(x = regions, crs = st_crs(countries))
regions$ID <- as.numeric(as.factor(regions$DEN_UTS))
plot(regions)

##--- Alps
alps <- read_sf(file.path(dataDir,"GISData/Output_layout/Italian_Alps.shp"))
alps <- st_transform(x = alps, crs = st_crs(countries))
plot(alps)

##---- Plot check
pdf(file = file.path(thisDir, "figures", paste0(modelName, "_ngs_map.pdf" )),
    width = 15, height = 12)
cols <- met.brewer(name="Isfahan1",n=10,type="continuous")
plot(st_geometry(st_intersection(studyArea,countries)),border = F)
plot(st_geometry(countries), col = "gray80", add = T, border = F)
plot(st_geometry(temp), col = "gray80",add=T,border = F)
plot(st_geometry(st_intersection(studyArea,countries)),add=T,col="gray60",border = F)
plot(transects,add=T,col="red")
plot(ngs[ngs$Sex == "M", ], add=T, cex = 1, col = cols[5], bg = adjustcolor(cols[5],0.3),pch=21)
plot(ngs[ngs$Sex == "F", ], add=T, cex = 1, col = cols[7], bg = adjustcolor(cols[7],0.3),pch=21)
legend( "bottomright", pch = 21, cex = 1.5, bty = "n", 
        legend = c("female","male"),
        col = cols[c(5,7)], 
        pt.bg = c(adjustcolor(cols[5],0.3),
                  adjustcolor(cols[7],0.3)))
graphics.off()


plot(st_geometry(countries), col = "gray80", border = F)
plot(st_geometry(temp), col = "gray80",add=T,border = F)
plot(st_geometry(st_intersection(studyArea,countries)),add=T,col="red",border = F)

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
plot(st_geometry(temp),col="gray60")
temp <- st_transform(temp, st_crs(studyArea))


## ------   2. SEARCH EFFORT DATA ------
##---- Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))

##---- Convert dates
transects$Date <- parse_date_time(transects$date, orders = c('ymd'))
transects$Year <- as.numeric(format(transects$Date,"%Y"))
transects$Month <- as.numeric(format(transects$Date,"%m"))

##---- Plot check
plot(transects, col = "red", add = T)



## ------   3. PRE-PROCESSED STUFF ------ 
load(file.path(thisDir,"Habitat.RData"))
load(file.path(thisDir,"Detectors.RData"))

# Check
# plot(detectors$grid)
# plot(habitat$grid)
# summary(habitat$grid)
# summary(detectors$grid)
# 
# plot(detectors$grid)
# plot(habitat$grid)

## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
## ------   1. RESCALE HABITAT & DETECTORS ------
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





## ------   2. DNA DATA ------
##---- All collected NGS samples
# allSamples <- read_sf(file.path(dataDir,"GISData/scats/merged_scats.shp"))
# plot(st_geometry(studyArea))
# plot(st_geometry(countries), col = "gray80",add=T)
# plot(allSamples, add = T,pch=3)

##---- Genotyped NGS samples
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

##---- Subset ngs data 
ngs <- ngs %>% sample_frac(0.5)
dim(ngs)
 
# ##---- Number of detections per individual
# numDetsPerId <- table(ngs$Genotype.ID)
# hist(numDetsPerId)
# 
# ##---- Number of individuals detected
# length(numDetsPerId)
# 
# ##---- Mean number of detections per ID
# mean(numDetsPerId)
# 
# ##---- Number of individuals with > 1 detection
# sum(numDetsPerId > 1)
# 
# ##---- mean number of individual detections by sex
# table(ngs$Sex, useNA = "always")
# apply(table(ngs$Genotype.ID, ngs$Sex, useNA = "always"),
#       2,
#       function(x)sum(x)/sum(x>0))
# 
# ##---- mean number of individual detections by status
# ngs$Status[ngs$Status == "pup 2016"] <- "pup"
# ngs$Status[ngs$Status == "pup 2017"] <- "pup"
# table(ngs$Status, useNA = "always")
# apply(table(ngs$Genotype.ID, ngs$Status, useNA = "always"),
#       2,
#       function(x)sum(x)/sum(x>0))
# 
# ##---- mean number of individual detections by sex and status
# apply(table(ngs$Genotype.ID, ngs$Status, ngs$Sex, useNA = "always"),
#       c(2,3),
#       function(x)sum(x)/sum(x>0))

##---- Turn ngs into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$CoordX, ngs$CoordY)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

## ------   3. DETECTION DATA ------
## ------     3.1. DETECTION MATRIX : y ------
##---- Calculate distance between detections and detectors
# closest <- nn2( st_coordinates(st_centroid(detectors$grid)),
#                 st_coordinates(ngs),
#                 k = 1,
#                 searchtype = "radius",
#                 radius = 10000)


# closest_ngs <- vector(mode = "list", length(ngs_ss))
# 
# for (i in 1:seq_along(ngs_ss)){
#   ngs_data <- ngs_ss[i]
#   closest <- nn2(coordinates(detectors$sub.sp),
#                   st_coordinates(ngs_data),
#                   k = 1,
#                   searchtype = "radius",
#                   radius = 10000)
#   closest_ngs[i] <- closest
# }
#   
closest <- nn2(coordinates(detectors$sub.sp),
                 st_coordinates(ngs),
                 k = 1,
                 searchtype = "radius",
                 radius = 10000)


##---- Assign each detection to a detector based on minimum distance
ngs$sub.detector <- c(closest$nn.idx)
ngs$detector <- detectors$sub.sp$main.cell.new.id[closest$nn.idx] 

# main.cell.id <- detectors$sub.sp$main.cell.new.id[closest$nn.idx]     
# ngs$detector <- unlist(
#   lapply(1:length(main.cell.id), function(x){
#     out <- NA
#     if(!is.na(main.cell.id[x]))out <- which(detectors$sp$main.cell.new.id %in% main.cell.id[x])
#     return(out)
#   }))
# all(ngs$detector==main.cell.id)

##---- Drop duplicated detections ate the same sub-detectors
ngs <- ngs[!duplicated(ngs[ ,c("sub.detector", "Genotype.ID")]), ] 
# ngs <- droplevels(ngs)

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


## ------     3.2. INDIVIDUAL COVARIATES ------
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




## ------     3.3. DATA AUGMENTATION ------


yCombined.aug <- MakeAugmentation( y = yCombined,
                                   M = 3000,
                                   replace.value = 0)

sex.aug <- MakeAugmentation( y = sex,
                             M = 3000,
                             replace.value = NA)

status.aug <- MakeAugmentation( y = status,
                                M = 3000,
                                replace.value = NA)



## -----------------------------------------------------------------------------
## ------ III. NIMBLE ------- 
## ------   1. MODEL ------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS  
  for(c in 1:n.habCovs){
    betaHab[c] ~ dnorm(0.0,0.01)
  }#c
  habIntensity[1:n.habWindows] <- exp(
    hab.covs[1:n.habWindows,1:n.habCovs] %*% betaHab[1:n.habCovs])
  
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
  
  for(s in 1:n.states){
    for(ss in 1:2){
      p0[s,ss] ~ dunif(0,0.5)
      sigma[s,ss] ~ dunif(0,12)
      logit(p0Traps[s,ss,1:n.detectors]) <- logit(p0[s,ss]) + 
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
                   "zone" = detectors$grid$`zone`,
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

nimParams <- c("N", "p0", "sigma", "psi",
               "betaDet", "betaHab", "theta", "rho",
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
                    "betaHab" = rep(0,nimConstants$n.habCovs),
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
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
for(c in 1:1){
  print(system.time(
    runMCMCbites( mcmc = Cmcmc,
                  bite.size = 500,
                  bite.number = 4,
                  path = file.path(thisDir, paste0("output/chain",c)))
  ))
}

##---- Collect multiple MCMC bites and chains
nimOutput_noZ <- collectMCMCbites( path = file.path(thisDir, "output"),
                                   burnin = 0,
                                   param.omit = c("s","z","sex","status"))

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput_noZ)
graphics.off()


##---- Process and save MCMC samples
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 0)
res <- ProcessCodaOutput(nimOutput)

##---- Save processed MCMC samples
save(res, file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))


