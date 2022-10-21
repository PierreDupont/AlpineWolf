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
# library(sf)
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
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
modelName = "AlpineWolf.SubSample.TransectsRepetitions"
thisDir <- file.path(analysisDir, modelName)
InDir <- file.path(thisDir, "input")


## ------   4. FIT MODEL -----
sim_names <- c(3,6)

for (sc in sim_names) {
  for (rp in 1:10) {
    for (it in 1:4) {
      ##---- Load input file
      thisFile <- paste0("AlpineWolf.SubSample.TransectsRepetitions_",sc,"_",rp,"_",it,"rep.RData")
      print(thisFile)
      load(file.path(InDir,thisFile))
      
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
      outFileName <- paste0()
      print(system.time(
        runMCMCbites( mcmc = Cmcmc,
                      bite.size = 50,
                      bite.number = 2,
                      path = file.path(thisDir,"output",thisFile))))
    }#sc 
  }#rp
}#it

# ##---- Collect multiple MCMC bites and chains
# nimOutput_noZ <- collectMCMCbites( path = file.path(thisDir, "output"),
#                                    burnin = 0,
#                                    param.omit = c("s","z","sex","status"))
# 
# ##---- Traceplots
# pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
# plot(nimOutput_noZ)
# graphics.off()
