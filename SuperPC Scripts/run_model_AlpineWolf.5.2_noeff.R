
# ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

# read the c parameter from command line
# Usage:
# Rscript run_model_0.3_RJMCM_x.R C_VALUE

args = commandArgs(trailingOnly = TRUE)
c_value = args[1]



## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)
library(R.utils)


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = "../Source", modifiedOnly = F)

modelName = "AlpineWolf.5.2_noeff"

##---- Set-up directories
inPath <- file.path("input")       ## Directory with input files as .RData
# procPath <- file.path("processed") ## Directory to store processed input files
outPath <- file.path("output")     ## Directory to store processed output files

# ##---- Load one input file at random
# file.list <- list.files(inPath)
# 
# if (length(file.list) == 0) {print("all models are processed"); exit()}
# 
# set <- sample(file.list, 1)
# load(paste(inPath, set, sep = "/"))
# print(set)
# 
# ##---- Move the selected input file to the 'procPath' directory (to avoid repeated runs).
# file.rename( from = file.path(inPath, set),
#              to = file.path(procPath, set))                                 
load( file.path(inPath, paste0(modelName, "_", c_value, ".RData")))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
# modelName = "AlpineWolf.SubSample.NGS"
# thisDir <- file.path(analysisDir, modelName)
# InDir <- file.path(thisDir, "input")

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
                       monitors = c("N","p0","sigma","psi","betaDet","betaHab","theta","rho","z","s","status","sex"),
                       thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc


##---- Run nimble MCMC in multiple bites
print(system.time(
  runMCMCbites( mcmc = Cmcmc,
                bite.size = 500,
                bite.number = 40,
                path = file.path(outPath, paste0(modelName, "_", c_value, ".RData")),
                save.rds = TRUE)))

