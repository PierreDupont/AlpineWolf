################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- AlpineWolf.SubSample.rep_100 --------------------- #####
##### ---------------------- DATA SUB-SAMPLING --------------------------- #####
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
modelName = "AlpineWolf.SubSample.rep_100"
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
## ------ III. SUB-SAMPLE DATA ------
ngs <- as.data.frame(ngs)


for(frac in sim_names){
  print(paste0("####   Sub-sampling NGS data with proportion ", frac*100, "%"))
  
  for (rep in 1:max.rep){
    
    ngs_sub <- ngs %>% sample_frac(frac)
    
    ##---- Drop duplicated detections ate the same sub-detectors
    ngs_sub <- ngs_sub[!duplicated(ngs_sub[ ,c("sub.detector", "Genotype.ID")]), ]
    ngs_sub <- droplevels(ngs_sub)
    
    ##---- Count individual detections per detector
    detMat_sub <- as.matrix(table(ngs_sub[ , c("Genotype.ID", "detector")]))
    
    ##---- Retrieve the number of detectors at which each individual is detected
    detNums_sub <- apply(detMat_sub, 1, function(x) sum(x>0))
    
    ##--- Set-up matrices to store individual detection frequencies and indices
    n.detected_sub <- dim(detMat_sub)[1]
    detIndices_sub <- matrix(-1, n.detected_sub, max(detNums_sub)*2)
    ySparse_sub <- matrix(-1, n.detected_sub, max(detNums_sub)*2)
    
    ##---- Fill in the matrix
    for(i in 1:n.detected_sub){
      ##-- Get where detections occur (detectors index)
      detIndices_sub[i,1:detNums_sub[i]] <- as.numeric(names(which(detMat_sub[i, ] > 0)))
      ##-- Get detection frequencies (always 1 in case of Bernoulli) 
      ySparse_sub[i,1:detNums_sub[i]] <- detMat_sub[i,which(detMat_sub[i, ] > 0)]
    }
    
    yCombined_sub <- cbind(detNums_sub, ySparse_sub, detIndices_sub)
    
    
    ## ------     3.2. INDIVIDUAL COVARIATES ------
    ##---- List detected indivdual names
    IDs_sub <- dimnames(yCombined_sub)[[1]]
    
    ##---- Set-up vectors to store individual covariates
    sex_sub <- status_sub <- pack_sub <- rep(NA, length(IDs_sub))
    for(i in 1:length(IDs_sub)){
      ##-- Sex 
      temp_sub <- unique(ngs_sub$Sex[ngs_sub$Genotype.ID %in% IDs_sub[i]])
      if(length(temp_sub)>1)warning(print(paste("ID:", IDs_sub[i], " sex:", temp_sub)))
      sex_sub[i] <- ifelse(length(temp_sub>1),temp_sub[1],temp_sub)
      
      ##- Social status
      temp_sub <- unique(ngs_sub$Status[ngs_sub$Genotype.ID %in% IDs_sub[i]])
      if(length(temp_sub)>1){
        #warning(print(paste("ID:", IDs_sub[i], " status:", temp_sub)))
        if(any(temp_sub %in% "na")){
          status_sub[i] <- temp_sub[!(temp_sub =="na")]
        }else{
          status_sub[i] <- NA
        }
      } else {
        status_sub[i] <- temp_sub
      }
      
      ##-- Pack membership
      temp_sub <- unique(ngs_sub$Pack[ngs_sub$Genotype.ID %in% IDs_sub[i]])
      if(any(temp_sub %in% "dispersal")){
        status_sub[i] <- "dispersal"
      }
    }
    
    ##----- Convert to numerical values 
    sex_sub[sex_sub == "F"] <- 0
    sex_sub[sex_sub == "M"] <- 1
    sex_sub[sex_sub == ""] <- NA
    sex_sub <- as.numeric(sex_sub)
    
    status_sub[status_sub == "alpha"] <- 1
    status_sub[status_sub == "pup"] <- 2
    status_sub[status_sub == "other"] <- 3
    status_sub[status_sub == "dispersal"] <- 3
    status_sub[status_sub == ""] <- NA
    status_sub <- as.numeric(status_sub)
    
    
    
    ## ------     3.3. DATA AUGMENTATION ------
    yCombined.aug_sub <- MakeAugmentation( y = yCombined_sub,
                                           M = 4000,
                                           replace.value = 0)
    
    sex.aug_sub <- MakeAugmentation( y = sex_sub,
                                     M = 4000,
                                     replace.value = NA)
    
    status.aug_sub <- MakeAugmentation( y = status_sub,
                                        M = 4000,
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
    
    
    
    ## ------   2. BUNDLE DATA ------
    ##---- Set model constants, data & parameter simulated values (==inits)
    nimData_sub <- list( y = yCombined.aug_sub,
                         z = c(rep(1,n.detected_sub),
                               rep(NA,dim(yCombined.aug_sub)[1]-n.detected_sub)),
                         sex = sex.aug_sub,
                         status = status.aug_sub,
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
    
    nimConstants_sub <- list( n.individuals = nrow(nimData_sub$y),
                              n.maxDets = ncol(nimData_sub$y),
                              n.habWindows = habitat$n.HabWindows,
                              n.detectors = detectors$n.detectors, 
                              n.habCovs = ncol(nimData_sub$hab.covs),
                              n.detCovs = ncol(nimData_sub$det.covs),
                              n.states = 3,
                              n.localIndicesMax = localObjects$numLocalIndicesMax,
                              y.max = dim(habitat$matrix)[1],
                              x.max = dim(habitat$matrix)[2])
    
    nimParams <- c("N","p0","sigma","psi","betaDet","betaHab","theta","rho")   ## Params with thin rate 1
    nimParams2 <- c("z","s","status","sex")                                    ## Params with thin rate 2
    
    
    
    ## ------   3. SAVE THE INPUT ------
    for(c in 1:4){
      s.init_sub <- matrix(NA, nimConstants_sub$n.individuals, 2)
      for(i in 1:n.detected_sub){
        if(detNums_sub[i] > 1){
          s.init_sub[i, ] <- colMeans(detectors$scaledCoords[detIndices_sub[i,1:detNums_sub[i]], ])
        } else {
          s.init_sub[i, ] <- detectors$scaledCoords[detIndices_sub[i,1:detNums_sub[i]], ] + rnorm(2,0,0.1)
        }
      }#i
      for(i in (n.detected_sub + 1):nimConstants_sub$n.individuals){
        s.init_sub[i, ] <- rbernppAC( n = 1,
                                      lowerCoords = nimData_sub$lowerHabCoords,
                                      upperCoords = nimData_sub$upperHabCoords,
                                      logIntensities = log(rep(1,habitat$n.HabWindows)),
                                      logSumIntensity = log(sum(rep(1,habitat$n.HabWindows))),
                                      habitatGrid = nimData_sub$habitatGrid,
                                      numGridRows = nrow(nimData_sub$habitatGrid),
                                      numGridCols = ncol(nimData_sub$habitatGrid))
      }#i
      
      
      sex.init_sub <- rbinom(n = nimConstants_sub$n.individuals, 1, prob = 0.5)
      sex.init_sub[!is.na(nimData_sub$sex)] <- NA
      
      status.init_sub <- rcat(n = nimConstants_sub$n.individuals, prob = c(0.5,0.45,0.05))
      status.init_sub[!is.na(nimData_sub$status)] <- NA
      
      z.init_sub <- rbinom(n = nimConstants_sub$n.individuals, 1, prob = 0.1)
      z.init_sub[!is.na(nimData_sub$z)] <- NA
      
      nimInits_sub <- list( "s" = s.init_sub,
                            "z" = z.init_sub,
                            "sex" = sex.init_sub,
                            "status" = status.init_sub,
                            "psi" = 0.5,
                            "rho" = 0.5,
                            "theta" = cbind(c(0.5,0.45,0.05),
                                            c(0.5,0.3,0.2)),
                            "betaDet" = rep(0,nimConstants_sub$n.detCovs),
                            "betaHab" = rep(0,nimConstants_sub$n.habCovs),
                            "p0" = cbind(c(0.1,0.1,0.05),
                                         c(0.1,0.1,0.05)),
                            "sigma" = cbind(c(1,1,2),
                                            c(1,1,2)))
      
      nimData <- nimData_sub
      nimConstants<- nimConstants_sub
      nimInits <- nimInits_sub
      
      save( modelCode,
            nimData,
            nimConstants,
            nimInits,
            nimParams,
            nimParams2,
            file = file.path(thisDir, "input", 
                             paste0(modelName, "_", round(frac*100),"_", rep, "_", c,  "rep.RData")))
    }#c
    
    print(rep)
  }#rep
}#frac


# for(frac in c(0.25,0.50,0.75)){
#   print(paste0("####   Renaming input data w/ proportion ", frac*100, "%"))
#   for(rep in 1:100){
#     for(c in 1:3){
#     file.rename(from = file.path(thisDir, "input", 
#                                  paste0(modelName, "_", frac,"_", rep, "_", c,  "rep.RData")),
#                 to = file.path(thisDir, "input", 
#                                paste0(modelName, "_", round(frac*100),"_", rep, "_", c,  "rep.RData")))
#     
#       }
#     print(rep)
#   }
# }



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
  resLista[[sc]] <- do.call(rbind,tempLista)
}#frac

##-- Give names to the list
names(resLista) <- sim_names*100

##-- Combine results into a dataframe and save as .csv
res <- do.call(rbind, resLista)
write.csv(res,
          file = file.path(thisDir, "results", paste0(modelName,".csv")))


## -----------------------------------------------------------------------------
## ------ II. PLOT SUBSAMPLING OUTPUTS ------
## ------   1. LOAD PROCESSED RESULTS -----
all_res <- read.csv(file.path(thisDir, "results.csv"))
def_res <- read.csv(file.path(thisDir, "res5_2.csv"))
def_res["scenario"] <- "default"

all_res <- all_res[ ,-1]
def_res <- def_res[ ,-1]

def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]



## ------   N   -----
g <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  facet_wrap(vars(stat), scales = "free_y") + 
  geom_point(data = def_res, size =2, alpha = 1.2, color = def_colors) +
  # scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="N - Wolves abundance", fill ="Search events", y = "N", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal") 

ggsave("transectsrepetition_N_comp.jpeg", dpi = 300)



## ------   p0   -----
p0 <- all_res[,13:18]
p0["stat"] <- all_res$stat
p0["scenario"] <- all_res$scenario
p0 <- melt(p0,id.vars=c("scenario","stat"))
# p01 <- p0 %>% 
#   filter(stat != "lci"| stat != "uci")

p0_52 <- def_res[,13:18]
p0_52["stat"] <- def_res$stat
p0_52["scenario"] <- def_res$scenario
p0_52 <- melt(p0_52,id.vars=c("scenario","stat"))

x_labs <- c("female alpha", "male alpha", "female pup",
            "male pup", "female other", "male other")

## ------   plots
g <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2, scale = "width") +
  geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = p0_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="p0 - baseline detectability", fill ="Search events", y = "p0", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_p0_comp.jpeg", dpi = 300)



## ------   sigma   -----
sigma <- all_res[,21:26] grep(pattern = "sigma",x = names(all_res))
sigma["stat"] <- all_res$stat
sigma["scenario"] <- all_res$scenario
sigma <- melt(sigma,id.vars=c("scenario","stat"))

sigma_52 <- def_res[,21:26]
sigma_52["stat"] <- def_res$stat
sigma_52["scenario"] <- def_res$scenario
sigma_52 <- melt(sigma_52,id.vars=c("scenario","stat"))

## ------   plots
g <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2, scale = "width") +
  geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = sigma_52, size =0.8, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="σ - scale parameter", fill ="Search events", y = "sigma", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_sigma_comp.jpeg", dpi = 300)



## ------   DetCov   -----
DetCov <- c("Transects Length", "Operator Experience", "Snow", "East/West", "Log Human Pop")

detcov <- all_res[,3:7]
detcov["stat"] <- all_res$stat
detcov["scenario"] <- all_res$scenario
detcov <- melt(detcov,id.vars=c("scenario","stat"))

detcov_52 <- def_res[,3:7]
detcov_52["stat"] <- def_res$stat
detcov_52["scenario"] <- def_res$scenario
detcov_52 <- melt(detcov_52,id.vars=c("scenario","stat"))

## ------   plots

# a <- ggplot(data = detcov_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   scale_x_discrete(labels= DetCov) +labs(title="Means", x ="sex-status", y = "Detection Covariates", fill = "Simulation")


g <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = detcov_52, size =0.8, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_detcov_comp.jpeg", dpi = 300)


## ------   HabCov   -----
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Human Pop","Wolf presence")

habcov <- all_res[,8:12]
habcov["stat"] <- all_res$stat
habcov["scenario"] <- all_res$scenario
habcov <- melt(habcov,id.vars=c("scenario","stat"))

habcov_52 <- def_res[,8:12]
habcov_52["stat"] <- def_res$stat
habcov_52["scenario"] <- def_res$scenario
habcov_52 <- melt(habcov_52,id.vars=c("scenario","stat"))




g <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = habcov_52, size =0.8, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_habcov_comp.jpeg", dpi = 300)

## ------   psi   -----
g <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = def_res,size =2, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y")  +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")


ggsave("transectsrepetition_psi_comp.jpeg", dpi = 300)


## ------   rho   -----
g <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = def_res, size =2, alpha = 1.2,color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y") +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_rho_comp.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   9. PIERRE'S HINTS -----
load(file.path(thisDir, "parm.df.RData"))
load(file.path(thisDir, "results.RData"))

## PLOTS N
{ 
  graphics.off()
  pdf(file = file.path(thisDir, "results_N.pdf"),
      width = 15, height = 8)
  
  species <- unique(parm.df$species)
  scenarios <- 1:8
  valNoise <- seq(-0.35,0.35, length.out = 8)
  par(mfrow = c(1,2))
  
  ## RELATIVE BIAS
  for(se in c("high", "low")){
    ylim <- 0.5
    plot(1,1, xlim = c(0.5,3.5), ylim = c(-ylim,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "RB(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(-ylim,ylim,length.out = 5),
         labels = seq(-ylim,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    
    p <- 1
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc  &
                                 myResults[[p]]$searchEffort == se, ]
        try(plot.violins2( dat.list = list(na.omit((temp$mean - temp$sim.N)/temp$sim.N)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  
  ## COEFFICIENT OF VARIATION
  for(se in c("high", "low")){
    ylim <- 0.15
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "CV(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se,  ]
        try(plot.violins2( dat.list = list(na.omit(temp$CV)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  ## COVERAGE
  for(se in c("high", "low")){
    ylim <- 1
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "95%CI coverage", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0.95, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se, ]
        
        coverage <- mean(temp$lci <= temp$sim.N & temp$uci >= temp$sim.N, na.rm =T)
        try(points( y = coverage,
                    x = sp + valNoise[sc],
                    pch = 21, cex = 2,
                    bg = colors[sc],
                    col = myCols[sc]), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 0.5, y = 0.35,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  graphics.off()
}




## -----------------------------------------------------------------------------
