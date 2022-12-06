###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------

library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(dplyr)
library(MASS) 
library(reshape2) 
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


## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
# sourceDirectory(path = "/archive/home/ofriard/projects/virginia-boiani/AlpineWolf/Source", modifiedOnly = F)

sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))

## -----------------------------------------------------------------------------
## ------ 5. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 

modelName = "AlpineWolf.SubSample.rep_100"
thisDir <- file.path(analysisDir, modelName)



## -----------------------------------------------------------------------------
## ------   6. PROCESS OUTPUTS -----

inDir <- file.path(thisDir,"input/")



resLista <-list()  #create a list which will contain all your results

sim_names <- c("25","50","75")  #Here the names of your simulation  (i.e. 25,50,75,100)
rep <- 4
iter <- 1   
     for (sc in sim_names) {
       tempLista <-list() #temporary list for each scenario
           
      
       for (rp in 1:rep) {
           for ( it in iter) {

    ##---- Load input file
    inputs <- load(file.path(inDir, paste0("AlpineWolf.SubSample.rep_100_", sc, "_", rp,"_", it,"rep.RData")))
    print(inputs)

  
    # Count sexes
    sex <- as.integer(nimData$sex)
    sex[sex == 0] <- "F"
    sex[sex == 1] <- "M"
    sex[sex == NA] <- ""
    ord <- seq(1:(length(sex)))
    df <- data.frame(ord,sex)
      
    sex_count <- df %>%
        count(sex)  %>%
        na.omit() 

  
    # Count Statuses
    status <- as.integer(nimData$status)
    status[status == 1] <- "RI"
    status[status == 2] <- "offspring"
    status[status == 3] <- "other"
    status[status == NA] <- ""
    ord <- seq(1:(length(status)))
    df <- data.frame(ord,status)
      
      
    status_count <- df %>%
      count(status) %>%
      na.omit()
    
    
    # merge all stats in one df
    colnames(status_count) <- c('variable','n')
    colnames(sex_count) <- c('variable','n')
    input_data <- rbind(sex_count,status_count)

    
    
    # Count id and detections
    det <- as.data.frame(nimData$y)
    det_sum <- sum(det$detNums_sub)
    det_max <- max(det$detNums_sub)
    id_count <- sum(det$detNums_sub > 0)
    
    id_det <- as.data.frame(rbind(det_sum, det_max, id_count))
    id_det <- tibble::rownames_to_column(id_det, "variable")
    id_det$variable <- c("NDetections","MaxDetections","IDs")
    id_det <- rename(id_det,c('n'='V1'))

    input_data <- rbind(input_data, id_det)
    
    
    # ## Add scenario name
    # input_data$scenario <- sim_names[sc]
    # input_data$iteration <- rep[rp]
    # # input_data$chain <- iter[it]
    
    # store df per each repetition in a temporary list
    tempLista[[rp]] <- input_data
    
    }#it
     }#rp
    
    # store all scenarios in one final list
    resLista[[sc]] <- tempLista
    
    }#sc

in25 <- do.call("rbind", resLista[[1]])
in25["scenario"] <- "25"
in50 <- do.call("rbind", resLista[[2]])
in50["scenario"] <- "50"
in75 <- do.call("rbind", resLista[[3]])
in75["scenario"] <- "75"

all_res <- rbind(in25, in50, in75)


write.csv(all_res,"prova.csv")



