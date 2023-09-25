###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------
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

## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)

modelName = "AlpineWolf.SubSample.Transects_prov_tr"
thisDir <- file.path(analysisDir, modelName)


## ------   6. PLOT RESULTS -----

OutDir <- file.path(thisDir, "output/")
## -----------------------------------------------------------------------------
## ------   6. PROCESS OUTPUTS -----


resLista <-list()  #create a list which will contain all your results


rpp <- 1
sim_names <- c("25","50","75")  #Here the names of your simulation  (i.e. 25,50,75,100)
rep_t <- c("6")

for(sc in 1:length(sim_names)){
  for(num in 1:length(rep_t)){
    
    
    tempLista <-list() #temporary list for each scenario
    
    for(rp in 1:rpp){  #This will be from 1:50
      outputs <- list.files(OutDir, paste0("AlpineWolf.SubSample.Transects_prov_tr", "_", sim_names[sc], "_", rp, "_", rep_t[num], "_"))
      print(outputs)

            nimOutput <- collectMCMCbites(path = file.path(OutDir, outputs),
                                     burnin = 10, #This will be 10 from now on
                                     param.omit = c("s","z","sex","status"))
      
      # obtain columns names
      params <- colnames(nimOutput[[1]]) 
      # Process mcmc
      tmp <- ProcessCodaOutput(nimOutput)
      # Select only mean, sd and rhat
      sd <- as.data.frame(unlist(tmp$sd))
      Rhat <- as.data.frame(unlist(tmp$Rhat))
      mean <- as.data.frame(unlist(tmp$mean))
      # bind them in df
      res <- as.data.frame(cbind(mean,sd,Rhat))
      res2 <- tibble::rownames_to_column(res, "params")
      # fix columns names
      res3 <- setnames(res2, old = c('params','unlist(tmp$mean)','unlist(tmp$sd)', 'unlist(tmp$Rhat)'), 
                       new = c('params','means','sd','rhat'))
      # transpose and again fix columns names 
      res4 <- as.matrix(res3)
      res5 <- t(res4)
      res5 <- as.data.frame(res5)
      colnames(res5) <- res5[1,]
      res5 <- res5[-1,] 
      res5 <- tibble::rownames_to_column(res5, "stat")
      
      # store df per each repetition in a temporary list
      tempLista[[rp]] <- res5
      
    }#rp
    
    # store all scenarios in one final list
    resLista[[sc]] <- tempLista
    
  }#rep_t
}#sc

#
res25 <- do.call("rbind",resLista[[1]])
res25["scenario"] <- "25"
res50 <- do.call("rbind",resLista[[2]])
res50["scenario"] <- "50"
res75 <- do.call("rbind",resLista[[3]])
res75["scenario"] <- "75"


all_res <- rbind(res25, res50, res75)
all_res["rep"] <- "6"
all_res["subsampling"] <- "Transects Spatial and Intensity"
res_tot_rhat <- rbind(filter(all_res, stat == "rhat"))

res_tot_rhat[,2:32] <- sapply(res_tot_rhat[,2:32],as.double)
res_tot_rhat<-res_tot_rhat[,-1] 



store <- list()

for (c in 1:ncol(res_tot_rhat)) {
  
  check <- vector()
  
  check <- res_tot_rhat[res_tot_rhat[,c] >= 1.1, ]
  store[[c]] <- check %>%
    group_by(scenario, subsampling, rep) %>%
    summarise(Count = n())
  
}
names(store) <- colnames(res_tot_rhat[,1:34]) 








