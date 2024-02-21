

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

modelName = "AlpineWolf_SC_prov.0.1"

##---- Set-up directories
thisDir <- file.path(project_dir, modelName)       ## Directory with current model


## -----------------------------------------------------------------------------
## ------ V. PROCESS NIMBLE OUTPUT ------
## ------   0. PROCESS MCMC CHAINS ------
##---- Collect multiple MCMC bites and chains
nimOutput <- collectMCMCbites( path = file.path(thisDir, "output"),
                               burnin = 10)

##---- Traceplots
pdf(file = file.path(thisDir, paste0(modelName, "_traceplots.pdf")))
plot(nimOutput)
graphics.off()


##---- Process and save MCMC samples
res <- ProcessCodaOutput(nimOutput$samples)
res_sxy <- ProcessCodaOutput(nimOutput$samples2)
save(res, res_sxy, 
     file = file.path(thisDir, paste0(modelName,"_mcmc.RData")))



