

# ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

## ------ SET REQUIRED WORKING DIRECTORIES ------

project_dir = "/archive/home/mboiani/AlpineWolf"
modelName = "AlpineWolf.5.4_RJMCMC_2324"
SourcePath <- "/archive/home/mboiani/AlpineWolf/Source"
## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)
library(R.utils)

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(SourcePath, modifiedOnly = F)


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
thisDir <- modelName
##---- Set-up directories
inPath <- file.path(modelName, "input")       ## Directory with input files as .RData
outPath <- file.path(modelName, "output")     ## Directory to store processed output files
thisDir <- file.path(modelName)       ## Directory with current model


## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
nimOutput_noZ <- collectMCMCbites( path = file.path(thisDir, "output"),
                                   burnin = 20,
                                   param.omit = c("s","z","sex","status"))

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput_noZ)
graphics.off()


##---- Process and save MCMC samples
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 20)
res <- ProcessCodaOutput(nimOutput)

##---- Save processed MCMC samples
save(res,nimOutput,
     nimOutput_noZ,
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))





