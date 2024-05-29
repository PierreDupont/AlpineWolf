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
modelName = "AlpineWolf.SubSample.Transects_prov_tr"
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
max.rep = 4                    ## Here the nmax number of repetitions(i.e. 100)
sim_names = c(0.25,0.50,0.75)  ## Here the names of your simulation  (i.e. 25,50,75,100)
rep_t =c(3,6)
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
# prov <- read_sf(file.path(dataDir,"GISData/Output_layout/northern_provinces.shp"))

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
trans_buf <- st_buffer(transects, dist = 500)

##---- Split NGS data into systematic and opportunistic samples
ngs_sys <- st_filter(ngs,trans_buf)
ngs_opp <- ngs[!ngs$Sample.ID %in% ngs_sys$Sample.ID, ]

##---- Assign each systematic smaple to a transect
ngs_sys$transect_ID <- unlist(lapply(st_intersects(ngs_sys, trans_buf), function(x)x[1]))
ngs_opp$transect_ID <- NA


##---- Loop over the different subsampling fractions 
for (frac in sim_names){
  for (rep in 1:max.rep) {
    
    ## ------   1. SUBSAMPLE TRANSECTS  -----   
    ##---- Sub-sample transects
    transectsToKeep <- trans %>%
      mutate(ID = row_number()) %>%  #Applying row_number function and obtaining the indexes
      group_by(trans$PROVINCIA)%>%  # Grouping for provinces
      slice_sample(prop = 0.75)    # subsampling given proportions of dataset                       
   
    pathsToKeep <- filter(transects, transect_i %in% transectsToKeep$ID_APP)
    
    ## ------   2. SUBSAMPLE TRANSECT REPETITIONS & NGS DATA -----   
    ##---- Subset repetitions  
    for (num in rep_t){
    
      trans3 <- pathsToKeep %>% 
      group_by(transect_i) %>%
      slice(sample(1:n(), pmin(3, n()))) 
      
    trans3_b <- filter(trans_buf, path_id %in% trans3$path_id)
    # trans3_b <- st_buffer(trans3, dist = 500)
    
    ##---- Extract ngs falling inside trans subsampled buffers
    ngs3 <- st_filter(ngs_sys, trans3_b) 
    ngs3 <- st_intersection(ngs3, trans3_b)
    ngs3 <- ngs3 %>% filter(Date == Date.1)
    ngs3 <- ngs3[,(1:20)]
    ngs3 <- rbind(ngs3,ngs_opp)
    
    
    ##---- Extract length and number of transects in each grid cell
    detectors3 <- detectors
    
    intersection3 <- st_intersection(detectors3$grid, trans3) %>%
      mutate(LEN = st_length(.),
             QI = .$Q.index) %>%
      st_drop_geometry() %>%
      group_by(id) %>%
      summarise(transect_L3 = sum(LEN),               ## Get total length searched in each detector grid cell
                transect_N3 = length(unique(Date)),   ## Get total number of visits in each detector grid cell
                transect_qi3 = mean(QI))              ## Get mean transects quality index for each detector grid cell
    
    ##---- Store in detector grid
    detectors3$grid <- detectors3$grid %>%
      left_join(intersection3, by = "id")
    detectors3$grid$transect_L3[is.na(detectors3$grid$transect_L3)] <- 0
    detectors3$grid$transect_N3[is.na(detectors3$grid$transect_N3)] <- 1
    detectors3$grid$transect_qi3[is.na(detectors3$grid$transect_qi3)] <- 0
    detectors3$grid$transect_L3 <- scale(detectors3$grid$transect_L3)
    detectors3$grid$mean_transect_L3 <- scale(detectors3$grid$transect_L3/detectors3$grid$transect_N3)
    
    ##---- Drop duplicated detections ate the same sub-detectors
    ngs3 <- as.data.frame(ngs3)
    ngs3 <- ngs3[!duplicated(ngs3[ ,c("sub.detector", "Genotype.ID")]), ]
    ngs3 <- droplevels(ngs3)
    
    ##---- Count individual detections per detector
    detMat3 <- as.matrix(table(ngs3[ , c("Genotype.ID", "detector")]))
    
    # ##---- Convert to binary detections (might keep later for binomial model)
    # detMat[detMat > 0] <- 1
    
    ##---- Retrieve the number of detectors at which each individual is detected
    detNums3 <- apply(detMat3, 1, function(x) sum(x>0))
    
    ##--- Set-up matrices to store individual detection frequencies and indices
    n.detected3 <- dim(detMat3)[1]
    detIndices3 <- matrix(-1, n.detected3, max(detNums3)*2)
    ySparse3 <- matrix(-1, n.detected3, max(detNums3)*2)
    
    ##---- Fill in the matrix
    for(i in 1:n.detected3){
      ##-- Get where detections occur (detectors index)
      detIndices3[i,1:detNums3[i]] <- as.numeric(names(which(detMat3[i, ] > 0)))
      ##-- Get detection frequencies (always 1 in case of Bernoulli) 
      ySparse3[i,1:detNums3[i]] <- detMat3[i,which(detMat3[i, ] > 0)]
    }
    
    ##---- Combine individual detection number, frequencies and indices
    yCombined3 <- cbind(detNums3, ySparse3, detIndices3)
    
    
    
    ## ------   2. INDIVIDUAL COVARIATES ------
    ##---- List detected indivdual names
    IDs3 <- dimnames(yCombined3)[[1]]
    
    
    ##---- Set-up vectors to store individual covariates
    sex3 <- status3 <- pack3 <- rep(NA, length(IDs3))
    for(i in 1:length(IDs3)){
      ##-- Sex 
      temp3 <- unique(ngs3$Sex[ngs3$Genotype.ID %in% IDs3[i]])
      if(length(temp3)>1)warning(print(paste("ID:", IDs3[i], " sex:", temp3)))
      sex3[i] <- ifelse(length(temp3>1),temp3[1],temp3)
      
      ##- Social status
      temp3 <- unique(ngs3$Status[ngs3$Genotype.ID %in% IDs3[i]])
      if(length(temp3)>1){
        #warning(print(paste("ID:", IDs25[i], " status:", temp25)))
        if(any(temp3 %in% "na")){
          status3[i] <- temp3[!(temp3 =="na")]
        }else{
          status3[i] <- NA
        }
      } else {
        status3[i] <- temp3
      }
      
      
      ##-- Pack membership
      temp3 <- unique(ngs3$Pack[ngs3$Genotype.ID %in% IDs3[i]])
      if(any(temp3 %in% "dispersal")){
        status3[i] <- "dispersal"
      }
    }
    
    
    ##----- Convert to numerical values 
    sex3[sex3 == "F"] <- 0
    sex3[sex3 == "M"] <- 1
    sex3[sex3 == ""] <- NA
    sex3 <- as.numeric(sex3)
    
    status3[status3 == "alpha"] <- 1
    status3[status3 == "pup"] <- 2
    status3[status3 == "other"] <- 3
    status3[status3 == "dispersal"] <- 3
    status3[status3 == ""] <- NA
    status3 <- as.numeric(status3)
    
    
    
    ## ------   3. DATA AUGMENTATION ------
    yCombined.aug3 <- MakeAugmentation( y = yCombined3,
                                        M = 4000,
                                        replace.value = 0)
    
    sex.aug3 <- MakeAugmentation( y = sex3,
                                  M = 4000,
                                  replace.value = NA)
    
    status.aug3 <- MakeAugmentation( y = status3,
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
    nimData3 <- list( y = yCombined.aug3,
                      z = c(rep(1,n.detected3),
                            rep(NA,dim(yCombined.aug3)[1]-n.detected3)),
                      sex = sex.aug3,
                      status = status.aug3,
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
                        "transect_L" = detectors3$grid$transect_L3,
                        "transect_qi" = detectors3$grid$transect_qi3,
                        "snow_fall" = detectors3$grid$`snow_fall`,
                        "zone" = detectors3$grid$`zone`,
                        "log_pop" = detectors3$grid$`log_pop`),
                      size = detectors3$size,
                      detCoords = detectors3$scaledCoords,
                      localTrapsIndices = localObjects$localIndices,
                      localTrapsNum = localObjects$numLocalIndices,
                      habitatGrid2 = habitat$matrix,
                      habitatGrid = localObjects$habitatGrid)
    
    
    nimConstants3 <- list( n.individuals = nrow(nimData3$y),
                           n.maxDets = ncol(nimData3$y),
                           n.habWindows = habitat$n.HabWindows,
                           n.detectors = detectors3$n.detectors, 
                           n.habCovs = ncol(nimData3$hab.covs),
                           n.detCovs = ncol(nimData3$det.covs),
                           n.states = 3,
                           n.localIndicesMax = localObjects$numLocalIndicesMax,
                           y.max = dim(habitat$matrix)[1],
                           x.max = dim(habitat$matrix)[2])
    
    
    
    nimParams <- c("N", "p0", "sigma", "psi",
                   "betaDet", "betaHab", "theta", "rho",
                   "z", "s", "status", "sex")
    
    
    
    ## ------   6. SAVE THE INPUT ------
    for(c in 1:4){
      s.init3 <- matrix(NA, nimConstants3$n.individuals, 2)
      for(i in 1:n.detected3){
        if(detNums3[i] > 1){
          s.init3[i, ] <- colMeans(detectors$scaledCoords[detIndices3[i,1:detNums3[i]], ])
        } else {
          s.init3[i, ] <- detectors$scaledCoords[detIndices3[i,1:detNums3[i]], ] + rnorm(2,0,0.1)
        }
      }#i
      for(i in (n.detected3 + 1):nimConstants3$n.individuals){
        s.init3[i, ] <- rbernppAC( n = 1,
                                   lowerCoords = nimData3$lowerHabCoords,
                                   upperCoords = nimData3$upperHabCoords,
                                   logIntensities = log(rep(1,habitat$n.HabWindows)),
                                   logSumIntensity = log(sum(rep(1,habitat$n.HabWindows))),
                                   habitatGrid = nimData3$habitatGrid,
                                   numGridRows = nrow(nimData3$habitatGrid),
                                   numGridCols = ncol(nimData3$habitatGrid))
      }#i
      
      sex.init3 <- rbinom(n = nimConstants3$n.individuals, 1, prob = 0.5)
      sex.init3[!is.na(nimData3$sex)] <- NA
      
      status.init3 <- rcat(n = nimConstants3$n.individuals, prob = c(0.5,0.45,0.05))
      status.init3[!is.na(nimData3$status)] <- NA
      
      z.init3 <- rbinom(n = nimConstants3$n.individuals, 1, prob = 0.1)
      z.init3[!is.na(nimData3$z)] <- NA
      
      nimInits3 <- list( "s" = s.init3,
                         "z" = z.init3,
                         "sex" = sex.init3,
                         "status" = status.init3,
                         "psi" = 0.5,
                         "rho" = 0.5,
                         "theta" = cbind(c(0.5,0.45,0.05),
                                         c(0.5,0.3,0.2)),
                         "betaDet" = rep(0,nimConstants3$n.detCovs),
                         "betaHab" = rep(0,nimConstants3$n.habCovs),
                         "p0" = cbind(c(0.1,0.1,0.05),
                                      c(0.1,0.1,0.05)),
                         "sigma" = cbind(c(1,1,2),
                                         c(1,1,2)))
      
      nimData <- nimData3
      nimConstants<- nimConstants3
      nimInits <- nimInits3
      
      save( modelCode,
            nimData,
            nimConstants,
            nimInits,
            nimParams,
            file = file.path(thisDir, "input",
                             paste0(modelName, "_", frac*100, "_", rep, "_", num,"_", c, ".RData")))
    }#c
    print(paste0("Dataset ",rep, " with ", num, " transect repetitions is done!"))
  }#rep
}#num
}#frac




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
  for(num in 1:length(rep_t)){
    
  tempLista <- list() ## Temporary list for each scenario
  
  for(rep in 1:max.rep){  
    
    ## List output files for scenario "num" and repetition "rep" only
    outPattern <- paste0(modelName, "_",  sim_names[frac], "_", rep, "_", rep_t[num], "_") 
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
      
      ## Store only the minimum in an object to summarize your results (e.g. N[num,rep] <- mean(nimOutput[ ,"N"])
      ## merge all stats in one df
      res <- as.data.frame(rbind(means,sd,CV,lci,uci))
      ## create columns for stats
      res2 <- tibble::rownames_to_column(res, "stat")
      params2 <- append(params, "stat", 0)
      colnames(res2) <- params2
      
      ## Add scenario name
      res2$scenario <- sim_names[frac]
      res2$repetition <- rep_t[num]
      
      ## store df per each repetition in a temporary list
      tempLista[[rep]] <- res2
    }
  }#rep
  ## Store all scenarios in one final list
  resLista[[frac]] <- do.call(rbind,tempLista)
}#sc
}#rep_t
##-- Give names to the list
names(resLista) <- sim_names

##-- Combine results into a dataframe and save as .csv
res <- do.call(rbind, resLista)
write.csv(res,
          file = file.path(thisDir, "results", paste0(modelName,".csv")))


## -----------------------------------------------------------------------------  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
