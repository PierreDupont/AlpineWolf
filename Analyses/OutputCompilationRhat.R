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
library(tidyr)
library(data.table)


## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")
modelName = "AlpineWolf.SubSample.Transects_prov_tr"
thisDir <- file.path(analysisDir, modelName)
OutDir <- file.path(thisDir, "output")


## ------ 4. SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)




## -----------------------------------------------------------------------------
## ------ 5. PROCESS OUTPUTS -----
sim_names <- c(25,50,75)  ## names of your scenarios  (i.e. 25,50,75,100)
rep_t <- c(3,6)           ## names of the transect repetitions (i.e. 3,6,all)
rpp <- 1                  ## number of replicates per scenario

resLista <- list()  ## list which will contain all your results
for(sc in 1:length(sim_names)){
  tmpLista1 <- list() ## temporary list for each scenario
  
  for(num in 1:length(rep_t)){
    tmpLista2 <- list() ## temporary list for each replicate
    
    for(rp in 1:rpp){  
      outputs <- list.files(OutDir,
                            paste0("AlpineWolf.SubSample.Transects_prov_tr", "_", sim_names[sc], "_", rp, "_", rep_t[num], "_"))
      print(outputs)
      
      ##-- Compile MCMC bites
      nimOutput <- collectMCMCbites(
        path = file.path(OutDir, outputs),
        burnin = 10,                            ## This will be 10 from now on
        param.omit = c("s","z","sex","status")) ## We could remove "betaHab" and "betaDet" here to make things faster

      ##-- Process MCMC output
      tmp <- ProcessCodaOutput(nimOutput)
      
      ##-- Select only mean, sd and rhat
      res <- rbind("mean" = unlist(tmp$mean),
                   "sd" = unlist(tmp$sd),
                   "Rhat" = unlist(tmp$Rhat))
      res <- tibble::rownames_to_column(as.data.frame(res), "stat")

      ##-- Save scenario info
      res$trsct_scenario <- sim_names[sc]
      res$trsct_repetitions <- rep_t[num]
      res$replicate <- rpp
        
      ##-- Store df for each scenario and each repetition 
      tmpLista1[[rp]] <- res
    }#rp
    
    ##-- Rbind all results for transect repetitions "num"
    tmpLista2[[num]] <- do.call(rbind,tmpLista1)
  }#num
  
  ##-- Rbind all results for scenario "sc"
  resLista[[sc]] <- do.call(rbind,tmpLista2)
  }#sc

##-- Combine all results together
all_res <- do.call(rbind, resLista)
all_res["subsampling"] <- "Transects Spatial and Intensity"

##-- Separate "mean", "sd" and "Rhat"
all_Rhat <- filter(all_res, stat == "Rhat")
all_mean <- filter(all_res, stat == "mean")
all_sd <- filter(all_res, stat == "sd")



## -----------------------------------------------------------------------------
## ------ 6. CHECK Rhat -----
##-- Identify parameter names
param_names <- names(all_Rhat)[2:32]

##-- Check R-hat for each parameter separately
store <- lapply(param_names, function(c){
  all_Rhat[all_Rhat[ ,c] >= 1.1, ] %>% 
    group_by(trsct_scenario,
             trsct_repetitions,
             subsampling) %>%
    summarise(Count = n())
})
names(store) <- param_names




