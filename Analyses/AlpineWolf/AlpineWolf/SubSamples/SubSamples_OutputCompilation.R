###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------ 1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ 2. IMPORT REQUIRED LIBRARIES ------
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
library(fs)
library(purrr)
library(ggplot2)


## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ 4. SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))




## -----------------------------------------------------------------------------
## ------ 5. SET ANALYSIS CHARACTERISTICS -----
modelName = "AlpineWolf.SubSample.Transects_prov"

thisDir <- file.path(analysisDir, modelName)

load(file.path(thisDir, "ExtractionHabitat.RData"))
load(file.path(thisDir, "Habitat.RData"))

OutDir <- file.path(thisDir, "output/")


## -----------------------------------------------------------------------------
## ------ 6. PROCESS OUTPUTS -----

## Create a list which will contain all your results
resLista <- list()  

## Define loop
rpp <- 2
sim_names <- c("25", "50","75")  
# rep_t <- c("3")

## Loop over simulation scenarios
for(sc in 1:length(sim_names)) {
  
  ## Create a temporary list for each scenario
  tempLista <- list() 
  
  # ptm <- proc.time()
  
  ## Loop over repetitions
  for(rp in 1:rpp) {  
     # for (num in 1:length(rep_t)) {

    ## List output files for scenario "sc" and repetition "rp" only
     outputs <- list.files(OutDir, paste0(modelName,"_", sim_names[sc], "_", rp, "_"))
    # outputs <- list.files(OutDir, paste0(modelName,"_25_1_"))
    
    print(outputs)
 
    ## Collect bites from the different chains for this percentage and this repetition only
    nimOutput <- collectMCMCbites(
      path = file.path(OutDir, outputs),
      burnin = 20) 
    ## Proces output
    res <- ProcessCodaOutput(nimOutput)
    
    ## Extract density
    resLista_temp <- list()
    
    for(s in 0:1){
      for(ss in 1:3){
        thisStatus <- (res$sims.list$z == 1) &
          (res$sims.list$sex == s) &
          (res$sims.list$status == ss)
        
        tmp <- GetDensity(
          sx = res$sims.list$s[ , ,1],
          sy = res$sims.list$s[ , ,2],
          z = thisStatus,
          IDmx = habitat.id,
          aliveStates = 1,
          returnPosteriorCells = F,
          regionID = regions.rgmx)
        
        resLista_temp <- cbind(
          resLista_temp,
          c(mean(tmp$PosteriorAllRegions),
            sd(tmp$PosteriorAllRegions),
            cv(tmp$PosteriorAllRegions)/100,
            as.numeric(quantile(tmp$PosteriorAllRegions,0.975)),
            as.numeric(quantile(tmp$PosteriorAllRegions,0.025))))
      }#ss
    }#s
    
    tmp <- GetDensity(
      sx = res$sims.list$s[ , ,1],
      sy = res$sims.list$s[ , ,2],
      z = res$sims.list$z,
      IDmx = habitat.id,
      aliveStates = 1,
      returnPosteriorCells = F,
      regionID = regions.rgmx)
    
    resLista_temp <- cbind(
      resLista_temp,
      c(mean(tmp$PosteriorAllRegions),
        sd(tmp$PosteriorAllRegions),
        cv(tmp$PosteriorAllRegions)/100,
        as.numeric(quantile(tmp$PosteriorAllRegions,0.975)),
        as.numeric(quantile(tmp$PosteriorAllRegions,0.025))))
    
    status <- c("Alpha","Pups","Others")
    colnames(resLista_temp) <- c(paste("N", "F", status, sep = "_"),
                            paste("N", "M", status, sep = "_"),
                            "N")
    
    ## store df per each repetition in a temporary list
    tempLista[[rp]] <-  cbind.data.frame(
      "stat" = c("mean","sd","cv","uci","lci"),
      resLista_temp)
  }#rp
  
  # totalTime <- proc.time()-ptm
  
  ## store all scenarios in one final list
  resLista[[sc]] <- tempLista
  # }#num
  }#sc


# Extracts lists with scenarios stored
res25 <- do.call("rbind", resLista[[1]])
res25["scenario"] <- "25"
res50 <- do.call("rbind", resLista[[2]])
res50["scenario"] <- "50"
res75 <- do.call("rbind", resLista[[3]])
res75["scenario"] <- "75"

# merge
all_res <- rbind(res25, res50, res75)
# add extra column with number of repetitions
# all_res["rep"] <- "6"
# save
# write.csv(all_res, "results_transectsubsample_prov_tr6.csv")










