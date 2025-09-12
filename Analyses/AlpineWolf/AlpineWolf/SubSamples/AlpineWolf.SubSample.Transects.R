################################################################################
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
library(purrr)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
# sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
# sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.Transects"
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

## SUBSAMPLING SPECIFICATION
max.rep = 100                       ## Here the nmax number of repetitions(i.e. 100)
sim_names = c(0.25,0.50,0.75)  ## Here the names of your simulation  (i.e. 25,50,75,100)

## SET DIRECTORIES
if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}
if(!dir.exists(file.path(thisDir, "results"))){dir.create(file.path(thisDir, "results"))}


## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. SEARCH EFFORT DATA ------
##---- Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))
trans <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/transects_wolfalps_2020-2021.shp"))

##---- Convert dates
transects$Date <- parse_date_time(transects$date, orders = c('ymd'))
transects$Year <- as.numeric(format(transects$Date,"%Y"))
transects$Month <- as.numeric(format(transects$Date,"%m"))



## ------   2. PRE-PROCESSED STUFF ------ 
load(file.path(thisDir,"Habitat.RData"))
load(file.path(thisDir,"Detectors.RData"))
detectors$grid <- detectors$grid[,-(3:6)]



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

##---- Turn ngs_sys into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$CoordX, ngs$CoordY)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)



## ------   3. DETECTION DATA ------
##---- Calculate distance between detections and detectors
closest <- nn2( coordinates(detectors$sub.sp),
                st_coordinates(ngs),
                k = 1,
                searchtype = "radius",
                radius = 10000)

##---- Assign each detection to a detector based on minimum distance
ngs$sub.detector <- c(closest$nn.idx)
ngs$detector <- detectors$sub.sp$main.cell.new.id[closest$nn.idx] 






## -----------------------------------------------------------------------------
## ------ III. SUBSAMPLE TRANSECTS ------
##---- Buffer all transects with 500m radius
trans_buf <- st_buffer(trans, dist = 500)

##---- Split NGS data into systematic and opportunistic samples
ngs_sys <- st_filter(ngs,trans_buf)
ngs_opp <- ngs[!ngs$Sample.ID %in% ngs_sys$Sample.ID, ]

##---- Assign each systematic smaple to a transect
ngs_sys$transect_ID <- unlist(lapply(st_intersects(ngs_sys, trans_buf), function(x)x[1]))
ngs_opp$transect_ID <- NA

##---- Loop over the different subsampling fractions 
for (frac in sim_names){
  for (rep in 1:max.rep) {
    ## ------   1. SUBSAMPLE TRANSECTS & NGS DATA -----   
    ##---- Sub-sample transects
    transectsToKeep <- which(rbinom(n = nrow(trans),
                                    size = 1,
                                    prob = frac) == 1)
    
    ##---- Remove samples associated w/ sub-sampled transects
    ngs_sub <- filter(ngs_sys, transect_ID %in% transectsToKeep)
    
    ##---- Keep opportunistic samples
    ngs_sub <- rbind(ngs_sub,ngs_opp)
    ngs_sub <- as.data.frame(ngs_sub)
    
    ##---- Extract length and number of transects in each grid cell
    ids <- Reduce(intersect, list(trans[transectsToKeep, ]$ID_APP, transects$transect_i))
    transects_sub <- subset(transects, transect_i  %in% ids)
    
    intersection <- st_intersection(detectors$grid, transects_sub) %>%
      mutate(LEN = st_length(.),
             QI = .$Q.index) %>%
      st_drop_geometry() %>%
      group_by(id) %>%
      summarise(transect_L = sum(LEN),               ## Get total length searched in each detector grid cell
                transect_N = length(unique(Date)),   ## Get total number of visits in each detector grid cell
                transect_qi = mean(QI))              ## Get mean transects quality index for each detector grid cell
    
    ##---- Store in detector grid
    detectors_sub <- detectors
    
    detectors_sub$grid <- detectors_sub$grid %>%
      left_join(intersection, by = "id")
    detectors_sub$grid$transect_L[is.na(detectors_sub$grid$transect_L)] <- 0
    detectors_sub$grid$transect_N[is.na(detectors_sub$grid$transect_N)] <- 1
    detectors_sub$grid$transect_qi[is.na(detectors_sub$grid$transect_qi)] <- 0
    detectors_sub$grid$transect_L <- scale(detectors_sub$grid$transect_L)
    detectors_sub$grid$mean_transect_L <- scale(detectors_sub$grid$transect_L/detectors_sub$grid$transect_N)
    
    
    ##---- Drop duplicated detections ate the same sub-detectors
    ngs_sub <- ngs_sub[!duplicated(ngs_sub[ ,c("sub.detector", "Genotype.ID")]), ]
    ngs_sub <- droplevels(ngs_sub)
    
    ##---- Count individual detections per detector
    detMat <- as.matrix(table(ngs_sub[ , c("Genotype.ID", "detector")]))
    
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
    
    
    
    ## ------   2. INDIVIDUAL COVARIATES ------
    ##---- List detected individual names
    IDs <- dimnames(yCombined)[[1]]
    
    ##---- Set-up vectors to store individual covariates
    sex <- status <- pack <- rep(NA, length(IDs))
    for(i in 1:length(IDs)){
      ##-- Sex 
      temp <- unique(ngs_sub$Sex[ngs_sub$Genotype.ID %in% IDs[i]])
      if(length(temp)>1)warning(print(paste("ID:", IDs[i], " sex:", temp)))
      sex[i] <- ifelse(length(temp>1),temp[1],temp)
      
      ##- Social status
      temp <- unique(ngs_sub$Status[ngs_sub$Genotype.ID %in% IDs[i]])
      if(length(temp)>1){
        #warning(print(paste("ID:", IDs[i], " status:", temp)))
        if(any(temp %in% "na")){
          status[i] <- temp[!(temp =="na")]
        }else{
          status[i] <- NA
        }
      } else {
        status[i] <- temp
      }
      
      ##-- Pack membership
      temp <- unique(ngs_sub$Pack[ngs_sub$Genotype.ID %in% IDs[i]])
      if(any(temp %in% "dispersal")){
        status[i] <- "dispersal"
      }
    }
    
    ##----- Convert to numerical values 
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
    
    
    
    ## ------   3. DATA AUGMENTATION ------
    yCombined.aug <- MakeAugmentation( y = yCombined,
                                       M = 4000,
                                       replace.value = 0)
    
    sex.aug <- MakeAugmentation( y = sex,
                                 M = 4000,
                                 replace.value = NA)
    
    status.aug <- MakeAugmentation( y = status,
                                    M = 4000,
                                    replace.value = NA)
    
    
    ## ------   4. NIMBLE MODEL ------
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
          sigma[s,ss] ~ dunif(0,20)
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
    
    
    
    
    ## ------   5. BUNDLE DATA ------
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
                       "transect_L" = detectors_sub$grid$transect_L,
                       "transect_qi" = detectors_sub$grid$transect_qi,
                       "snow_fall" = detectors_sub$grid$`snow_fall`,
                       "zone" = detectors_sub$grid$`zone`,
                       "log_pop" = detectors_sub$grid$`log_pop`),
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
    
    
    
    ## ------   6. SAVE THE INPUT ------
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
      
      nimData <- nimData
      nimConstants<- nimConstants
      nimInits <- nimInits
      
      save( modelCode,
            nimData,
            nimConstants,
            nimInits,
            nimParams,
            file = file.path(thisDir, "input",
                             paste0(modelName, "_", frac*100, "_", rep, "_", c, ".RData")))
    }#c
    print(paste0("Dataset ",rep, " with ", frac*100, "% transect sub-sampling is done!"))
  }#rep
}#frac




## -----------------------------------------------------------------------------




## ------ IV. FIT MODEL -----
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
system.time(
  runMCMCbites( mcmc = Cmcmc,
                bite.size = 500,
                bite.number = 2,
                path = file.path(thisDir, paste0("output/chain",c)))
)


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



## -----------------------------------------------------------------------------
## ------ V. PROCESS OUTPUTS -----
##-- create a list which will contain all your results
resLista <- list()  
for(frac in 1:length(sim_names)){
  
  tempLista <- list() ## Temporary list for each scenario
  
  for(rep in 1:max.rep){  
    
    ## List output files for scenario "frac" and repetition "rep" only
    outPattern <- paste0(modelName, "_",  sim_names[frac]*100, "_", rep, "_") 
    outputs <- list.files(outDir,  pattern = outPattern)
    
    ## continue only if there are some outputs matching the pattern
    if(length(outputs) > 0){
      
      print(outputs)
      
      ## Check to make sure all files listed in "outputs" are bigger than 0 KB 
      if (any( file.size(file.path(outDir, outputs)) <= 0)) stop(paste0("One of ",outPattern," files is 0 KB!"))
      
      ## Collect bites from the different chains 
      nimOutput <- collectMCMCbites( path = file.path(outDir, outputs),
                                     burnin = 0,
                                     param.omit = c("s","z","sex","status"))
      
      ## Continue processing for this output inside the loop (e.g. extract mean values and CI for N, sigma etc...)
      ## Attach chains
      nimOutput <- do.call(rbind, nimOutput$samples)
      ## Obtain columns names
      params <- colnames(nimOutput)
      ## turn nimOutput as df 
      nimOutput <- as.data.frame(nimOutput)
      
      ## run functions to get stats
      means <- col_mean(nimOutput, params)
      sd <- col_sd(nimOutput,params)
      CV <- col_cv(nimOutput,params)
      uci <- col_uci(nimOutput,params)
      lci <- col_lci(nimOutput,params)
      
      ## Store only the minimum in an object to summarize your results (e.g. N[frac,rep] <- mean(nimOutput[ ,"N"])
      ## merge all stats in one df
      res <- as.data.frame(rbind(means,sd,CV,lci,uci))
      ## create columns for stats
      res2 <- tibble::rownames_to_column(res, "stat")
      params2 <- append(params, "stat", 0)
      colnames(res2) <- params2
      
      ## Add scenario name
      res2$scenario <- sim_names[frac]*100
      
      ## store df per each repetition in a temporary list
      tempLista[[rep]] <- res2
    }
  }#rep
  ## Store all scenarios in one final list
  resLista[[frac]] <- do.call(rbind,tempLista)
}#frac

##-- Give names to the list
names(resLista) <- sim_names*100

##-- Combine results into a dataframe and save as .csv
res <- do.call(rbind, resLista)
write.csv(res,
          file = file.path(thisDir, "results", paste0(modelName,".csv")))


## -----------------------------------------------------------------------------