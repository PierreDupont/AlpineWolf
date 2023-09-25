## -------------------------------------------------------------------------- ##
## ----------------------- ALPINE WOLF SCR ---------------------------------- ##
## ---------------------- CONVERSION FACTOR---------------------------------- ##
## -------------------------------------------------------------------------- ##
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
library(MetBrewer)
library(data.table)
library(ggplot2)

## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.5.2_RJMCMC"
thisDir <- file.path(analysisDir, modelName)

##---- Load processed MCMC samples
load(file = file.path(thisDir, paste0(modelName, "_densities.RData")))




## ------ 1. CALCULATE CONVERSION FACTOR FOR EAST/WEST -----
posterior_N <- WA_regions$PosteriorRegions
total <- colSums(posterior_N)
posterior_N <- rbind(posterior_N,total)
rownames(posterior_N) <- c("East","West","Total")

posterior_RI <- (WA_status[["F"]][["Alpha"]]$PosteriorRegions + WA_status[["M"]][["Alpha"]]$PosteriorRegions)/2
total <- colSums(posterior_RI)
posterior_RI <- rbind(posterior_RI,total)
rownames(posterior_RI) <- c("East","West","Total")

posterior_CF <- posterior_N/posterior_RI
apply(posterior_CF,1,mean)
apply(posterior_CF,1,function(q)quantile(q,c(0.025,0.5,0.975)))


posterior_CF2 <- as.data.frame(t(posterior_CF)) 


## ------ 2. VIOLIN PLOTS (all iterations) -----
# pdf(file.path(analysisDir, modelName, "Coefficient_Plot.pdf"),
#     width = 12, height = 10)
post_CF <- pivot_longer(posterior_CF2,c("East","West","Total"))

ggplot(post_CF, aes(value, name)) +
  geom_violin( draw_quantiles = c(0.025, 0.5, 0.975),
               fill = "darkgreen",
               color = "white")

# dev.off()

## -----------------------------------------------------------------------------