

# ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

# read the c parameter from command line
# Usage:
# Rscript run_model_SC_prov_0.1.R

args = commandArgs(trailingOnly = TRUE)
c_value = args[1]



## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)
library(R.utils)

project_dir = "/archive/home/ofriard/projects/virginia-boiani/SC_model"

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = "Source", modifiedOnly = F)

modelName = "AlpineWolf.5.2_noeff"

##---- Set-up directories
thisDir <- file.path(project_dir, modelName)       ## Directory with current model


## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
nimOutput_noZ <- collectMCMCbites( path = file.path(thisDir, "output"),
                                   burnin = 10,
                                   param.omit = c("s","z","sex","status"))

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput_noZ)
graphics.off()


##---- Process and save MCMC samples
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 10)
res <- ProcessCodaOutput(nimOutput)

##---- Save processed MCMC samples
save(res,nimOutput,nimOutput_noZ,
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))





